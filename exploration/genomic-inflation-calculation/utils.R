library(data.table)
library(tidyverse)
library(beepr)
library(ggpubr)
source('themes/theme_main.R')

generate_gi_lambda_data <- function(path, fig_out_path, alert_when_done = FALSE) {
      #' Generate the genomic inflation data for given compounds in the path as
      #' well as the supplementary figures for Manhattan and QQplots.
      #'
      #' @param path path to the folder containing the *_std_pvals.csv files
      #' @param fig_out_path path where Manhattan and QQ-plot figures should be stored
      #' @param alert_when_done if TRUE, an alert will sound when the code finishes running
      #' @return a data frame with the following columns:
      #'        - CompoundName, the name of the compound
      #'        - Lambda, the genomic inflation value

  # gather all the pvalue files in the path
  list_of_files <- list.files(path, pattern = "*std_pvals.csv", full.names = TRUE, recursive = TRUE)
  
  # create fig_out_path if does not exist
  if( !dir.exists(fig_out_path) ){
    dir.create(fig_out_path)
  }
  
  # create a final dataframe
  lambdas.df <- data.frame(matrix(ncol = 2, nrow = 0))

  for (file in list_of_files) {
    
    # evaluate the compound name from the URI
    compound_name <- sub(
      "_std_pvals.csv", "",
      sub(
        ".*/", "",
        sub(path, "", file)
      )
    )

    # read pvals, and SNP data from the file
    fs <- fread(file)
    pvals <- fs$pval
    snp <- fs$SNP
    
    # calculate lambda value for the compound
    chisq <- qchisq(pvals, 1)
    lmbd.val <- round(median(chisq) / qchisq(0.5, 1),3)
    
    # generate GWAS data used for Manhattan and qq plots
    gwas_data <- data.frame(SNP=snp,pvals=pvals) %>%
      separate(SNP,into = c("chr","BP","Tool","Allele"), sep = "_", remove = F) %>%
      mutate(CHR=str_sub(chr,2,-1)) %>%
      mutate(CHR=as.integer(CHR), BP=as.integer(BP)) %>%
      mutate(CHR=replace(CHR,CHR == 0, 18))
    
    # Compute chromosome size
    don <- gwas_data %>% group_by(CHR) %>%
      summarize(chr_len=max(BP)) %>%
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>% 
      left_join(gwas_data,., by=c("CHR"="CHR")) %>%
      arrange(CHR,BP) %>% 
      mutate(BPcum=BP+tot)
    
    # calculate the axis data
    axisdf <- 
      don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    
    # calculate the Bonferroni correction threshold for horizontal line in 
    # Manhattan plot
    bonf_thresh <- -log10(0.05 / nrow(don))
    
    # generate the Manhattan plot
    manhatt.plot <- ggplot(don, aes(x=BPcum, y=-log10(pvals))) +
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("#A4A4A4", "#2E2E2E"), 22 )) +
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0,0), limits = c(0,max(-log10(pvals)+2))) +
      GLOBAL_THEME + 
      labs(y=expression(-log[10](p)), x="Chromosome") + 
      theme(
        legend.position = "none"
      ) + 
      geom_hline(yintercept = bonf_thresh, colour = "red")
    
    
    # generate the QQ plot
    o <- -log10(sort(gwas_data$pvals, decreasing=FALSE))
    e <- -log10(ppoints(length(gwas_data$pvals)))
    annot_x <- min(e) + 0.5
    annot_y <- max(o) - 1
    qq.plot <- data.frame(
      o = o,
      e = e
    ) %>% ggplot(aes(x=e,y=o)) + geom_point() + 
    scale_x_continuous(limits=c(0,max(e))) + 
    scale_y_continuous(limits=c(0, max(o))) +
    GLOBAL_THEME + labs(
      x=expression(Expected~~-log[10](italic(p))),
      y=expression(Observed~~-log[10](italic(p)))
    ) + geom_abline(slope = 1, colour = "red") +
      annotate("text", x = annot_x, y = annot_y, size = 11, label = expr(paste(lambda,'=',!!lmbd.val)))

    # This above code returns following warning because the 'label' in annotate function does not evaluate
    # to a list, but it evaluates to a language type variable. FYI:
    #
    # Warning in is.na(x) :
    # is.na() applied to non-(list or vector) of type 'language'

    # annotate("text",x=min(e)+0.5,y=max(o)-1, size=11,label= expr(paste(lambda,'=',!!lmbd.val)))
    
    final_plot <- ggarrange(manhatt.plot, qq.plot, nrow=1, ncol=2)
    plot_name <- 
      paste0(fig_out_path,compound_name,".png")
    ggsave(final_plot, filename = plot_name, bg="white", width=15.5, height=5.5)
    
    # populate the final lambdas.df data frame
    lambdas.df <- rbind(
      lambdas.df,
      c(compound_name, as.numeric(lmbd.val))

    )
  }
  colnames(lambdas.df) <- c("CompoundName", "Lambda")

  # fix some of the compound names so that it matches the GCMS dataset
  lambdas.df[grep("Geranylacetone", lambdas.df$CompoundName), 'CompoundName'] <- "(E)-Geranylacetone"
  lambdas.df[grep("furanone", lambdas.df$CompoundName), 'CompoundName'] <- "5-ethyl-2(5H)-furanone"
  lambdas.df$CompoundName <- gsub("_", " ", lambdas.df$CompoundName)

  # make the lambda column numeric
  lambdas.df$Lambda <- as.numeric(lambdas.df$Lambda)

  if (alert_when_done) beep(sound = 3) # if asked to alert when done, do so

  return(lambdas.df)
}

generate_wilcoxon_data <- function(gcms_data_file) {
      #' Generatethe Shapiro-Wilk test W-statistic and p-value data for given compounds in GCMS data file.
      #' This function returns the W statistic and p-values for normality test for the given GCMS compounds
      #'
      #' @param gcms_data_file the path to the GCMS data file
      #' @return a data frame with the following columns:
      #'        - CompoundName, the name of the compound
      #'        - W, the W-statistic
      #'        - p, the p-value

  gcms_pheno_tbl <- read_excel(gcms_data_file)
  gcms_pheno_tbl.noaid <- gcms_pheno_tbl[, 2:ncol(gcms_pheno_tbl)]

  # perform the normality test
  test_results <- lapply(gcms_pheno_tbl.noaid, shapiro.test)

  final.df <- data.frame(matrix(ncol = 3, nrow = 0))
  for (compound_name in colnames(gcms_pheno_tbl.noaid)) {
    w <- test_results[[compound_name]]$statistic
    p <- test_results[[compound_name]]$p.value

    final.df <- rbind(
      final.df,
      c(compound_name, w, p)
    )
  }

  row.names(final.df) <- NULL
  colnames(final.df) <- c("CompoundName", "W", "p")
  return(final.df)
}
