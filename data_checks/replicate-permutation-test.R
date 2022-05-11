# Title     : Replicate Permutation Test
# Objective : This script compares a random distribution of samples with the
#             biological replicates to see if the biological replicates are
#             significantly correlated
# Created by: tayabsoomro
# Created on: 2021-04-10

library(dplyr)
library(stringr)
library(combinat)
library(ggplot2)
library(ggpubr)


# nursery ids that are biological replicates
s1 <- c(6166,1008,1196,2249,9304,9132,8087)
s2 <- c(6167,1009,1197,2250,9305,9134,8089)

random_perm_check <- function(pheno_tbl,out_file) {
  # get the indexes of all the rows where the apple id is any of the replicates
  idx_rem1 <- idx_rem2 <- NULL
  for (i in seq_along(s1)) {
    idx1 <- which(rownames(pheno_tbl) == s1[i])
    idx2 <- which(rownames(pheno_tbl) == s2[i])
    
    if(length(idx1)) {
      idx_rem1[i] <- idx1
    }
    
    if(length(idx2)) {
      idx_rem2[i] <- idx2
    }
    
  }
  
  # The nursery ids 1196, 1197 and 1009 were not standardized and therefore they
  # are not included in the standardized version of the phenotype table
  
  # remove the rows which belong to replicates
  pheno_tbl <- pheno_tbl[-na.omit(idx_rem1),]
  pheno_tbl <- pheno_tbl[-na.omit(idx_rem2),]
  
  common_compounds <- est <- pval <- NULL
  lens <- NULL
  for(i in 1:10000){
    # sample two at random from phenotype table
    smpl <- sample(nrow(pheno_tbl),2)
    
    # compound abundance values for both samples
    s1_dat <- as.numeric( pheno_tbl[ smpl[1], ] )
    s2_dat <- as.numeric( pheno_tbl[ smpl[2], ] )
    
    # get the common compounds
    common_compounds[i] <- sum(s1_dat != 0 & s1_dat != 0)
    
    # exclude the compounds that are not present in both samples
    idx_rem <- which(s1_dat == 0 & s2_dat == 0)
    s1_dat <- s1_dat[-idx_rem]
    s2_dat <- s2_dat[-idx_rem]
    
    # log transformation of the abundance values
    s1_dat_log <- log10(s1_dat)
    s2_dat_log <- log10(s2_dat)
    
    s1_dat_log[which(is.infinite(s1_dat_log))] <- NA
    s2_dat_log[which(is.infinite(s2_dat_log))] <- NA
    
    if( sum(!is.na(s1_dat_log)) > 3 & sum(!is.na(s2_dat_log)) > 3 ) {
      test <- cor.test(s1_dat_log,s2_dat_log)
      est[i] <- test$estimate
      pval[i] <- test$p.value
    }
    
  }

  # generate the histograms and save it to out_file
  png(out_file)
  par(mfrow=c(2,2))
  hist(
    common_compounds, 
    main="Distribution of \ncommon compounds",
    xlab="Number of samples",
    ylab="Number of compounds",
    breaks = 100
  )
  hist(
    est^2,
    main="Distribution of \nR² value",
    xlab="R²",
    ylab="Number of samples"
  )
  hist(
    pval,
    main="Distribution of \np-value",
    xlab="p-value",
    ylab="Number of samples",
    breaks = 100
  )
  dev.off()
}

# generate the permutation correlation graphs for random samples
pheno_tbl_hex <- 
  read.table('data/standardized/gcms_phenotype_table_std_w_hexanone.tsv')

random_perm_check(pheno_tbl_hex,'figs/standard-checks/perm_test_hex.png')

pheno_tbl_benz <- 
  read.table('data/standardized/gcms_phenotype_table_std_w_benzaldehyde.tsv')

random_perm_check(pheno_tbl_hex,'figs/standard-checks/perm_test_benz.png')

pheno_tbl_ethyl <- 
  read.table('data/standardized/gcms_phenotype_table_std_w_ethylacetate.tsv')

random_perm_check(pheno_tbl_hex,'figs/standard-checks/perm_test_ethyl.png')


# get the correlation graph of replicates
replicate_check <- function(pheno_tbl,standard_name,out_file) {
  plots <- list()
  for(i in seq_along(s1)){
    row1 <- which(rownames(pheno_tbl) == s1[i])
    row2 <- which(rownames(pheno_tbl) == s2[i])
    
    # get the numeric compound abundances
    dat1 <- as.numeric(pheno_tbl[row1,])
    dat2 <- as.numeric(pheno_tbl[row2,])
    
    print(max(dat2))
    
    # remove compounds that are not present in any of the samples
    idx_rem <- which(dat1 == 0 & dat2 == 0)
    dat1[-idx_rem]
    dat2[-idx_rem]
    
    # log transform the data
    dat1_log <- log10(dat1)
    dat2_log <- log10(dat2)
    
    
    plots[[i]] <- ggplot(data.frame(dat1_log,dat2_log),aes(x = dat1_log, y = dat2_log)) + 
      geom_point() + 
      xlab(paste0("Nursery ID: ", s1[i])) + 
      ylab(paste0("Nursery ID: ", s2[i])) + 
      theme_classic() +
      stat_cor(
        method = "spearman", 
        aes(label = paste(..r.label.., ..rr.label..,..p.label.., sep = "~`,`~"))
      )
    
  }
  
  fig <- annotate_figure(
    ggarrange(plotlist=plots, nrow = 3, ncol = 3),
    top=text_grob(
      paste0("Correlation of ", standard_name, "\nstandardized replicates"), 
      face = "bold", 
      size = 16
    )
  )
  
  png(out_file, width = 1000, height = 700)
  print(fig) 
  dev.off()
}

replicate_check(
  pheno_tbl_hex,
  "2-hexanone-1,1,1,3.3-d5",
  "figs/standard-checks/rep_test_hex.png"
)
replicate_check(
  pheno_tbl_benz,
  "Benzaldehyde-d6",
  "figs/standard-checks/rep_test_benz.png"
)
replicate_check(
  pheno_tbl_ethyl,
  "Ethyl acetate-d8",
  "figs/standard-checks/rep_test_ethyl.png"
)
