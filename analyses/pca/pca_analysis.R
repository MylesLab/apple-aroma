# Title     : PCA plots
# Objective : This script generates the PCA bi-plots as well as some of the other
#             plots made from PCA data

############################
## IMPORTS & DATA LOADING ##
############################

library(readxl)
library(ggpubr)
library(dplyr)
library(viridis)

source('analyses/pca/utils.R')
source('themes/theme_main.R')

# loading the GCMS phenotype table
gcms_pheno_tbl <-
  read_excel('data/processed/Supplementary_Data.xlsx', sheet = 'GCMS Data')
dim(gcms_pheno_tbl)
# [1] 515 107

# loading the ABC phenotype table
abc_pheno_tbl <- read_excel('data/processed/Supplementary_Data.xlsx', sheet = 'Apple Traits Data')
dim(abc_pheno_tbl)
# [1] 515 44

# loading the analysis date data
date_tbl <- openxlsx::read.xlsx('data/processed/harvest_date_analysis_date_pivot.xlsx')
date_tbl <- date_tbl[, c("AppleID", "AnalysisDateJulian", "HarvestDate")]
date_tbl$AppleID <- as.numeric(date_tbl$AppleID)

# merge the GCMS and ABC phenotype data
pheno_tbl <- left_join(gcms_pheno_tbl, abc_pheno_tbl, by = c("appleid" = "apple_id"))
dim(pheno_tbl)
# [1] 515 150


# add harvest date data to pheno_tbl
pheno_tbl <- inner_join(pheno_tbl, date_tbl, by = c("appleid" = "AppleID"))
dim(pheno_tbl)
# [1] 515 152

# get the column indices for aroma compounds
aroma_columns_idx <- seq_along(gcms_pheno_tbl)[-1]

# Gather the PCA data
pca_data <- pheno_tbl[, aroma_columns_idx]
dim(pca_data)
# [1] 515 106

##################
## PCA ANALYSIS ##
##################

# generate the pca data frame
pca_data.pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)
summary(pca_data.pca)

# calculate the ubiquity of volatiles for each sample
VolatileUbiquityBySample <- rowSums(pca_data != 0)

# calculate the abundance of volatiles for each sample
TotalVolatileAbundanceBySample <- rowSums(pca_data)


# calculate the variance of PCs
vars_transformed <- apply(pca_data.pca$x, 2, var)
pov <- (vars_transformed / sum(vars_transformed)) * 100

pcs.df <- data.frame(
  PC1 = pca_data.pca$x[, 1],
  PC2 = pca_data.pca$x[, 2],
  AppleId = pheno_tbl$appleid,
  Use = pheno_tbl$use,
  World = pheno_tbl$world,
  HarvestDate = pheno_tbl$date_jul_17_harv,
  VolatileUbiquityBySample = VolatileUbiquityBySample,
  TotalVolatileAbundanceBySample = TotalVolatileAbundanceBySample
)

lim <- c(
  min(pcs.df$PC1), max(pcs.df$PC1),
  min(pcs.df$PC2), max(pcs.df$PC2)
)
limits <- c(min(lim), max(lim))

########################
## GENERATE PCA PLOTS ##
########################

# pca plot for harvest date change
fig_3a.plot <- generate_pca_biplot(
  pcs.df,
  choices=c("PC1","PC2"),
  color_phenotype = c("HarvestDate", "Harvest Date\n (julian days)"),
  limits = limits,
  proportion_of_variance = pov
)
ggsave(
  "analyses/pca/figures/fig_3a.png",
  fig_3a.plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px"
)

# PCA bi-plot for powerpoint presentation
fig_3a.plot.ppt <- generate_pca_biplot(
  pcs.df,
  choices=c("PC1","PC2"),
  color_phenotype = c("HarvestDate", "Harvest Date\n (julian days)"),
  limits = limits,
  proportion_of_variance = pov,
  options = list(
    dot_size = 4,
    axis_title_size = 18,
    axis_text_size = 14,
    legend_title_size = 18,
    legend_text_size = 14
  )
)
ggsave(
  "analyses/pca/figures/Powerpoint_fig_3a.png",
  fig_3a.plot.ppt,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px"
)

# Figure 3B
fit.harv_dt_vs_pc1 <- lm(pcs.df$PC1 ~ pcs.df$HarvestDate)
summary(fit.harv_dt_vs_pc1)
fig_3b.plot <- ggplot(pcs.df, aes(x = HarvestDate, y = PC1)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_line(
    aes(x = HarvestDate, y = fit.harv_dt_vs_pc1$fitted.values),  # predicted data
    color = 'black') +
  geom_text(size=8,aes(label = sprintf("R² = %f", summary(fit.harv_dt_vs_pc1)$r.squared), x = 280, y = -15)) +
  geom_text(size=8,aes(label = sprintf("p-value = %s", format(summary(fit.harv_dt_vs_pc1)$coefficients[2, 4], scientific = TRUE)), x = 280, y = -17)) +
  GLOBAL_THEME +
  theme(
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24)
  ) +
  xlab("Harvest Date (julian days)") +
  ylab("PC1")
ggsave(
  "analyses/pca/figures/fig_3b.png",
  plot = fig_3b.plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px"
)

# Figure 3C
fit.harv_dt_vs_num_vols <- lm(pcs.df$VolatileUbiquityBySample ~ pcs.df$HarvestDate)
summary(fit.harv_dt_vs_num_vols)
fig_3c.plot <- ggplot(pcs.df, aes(x = HarvestDate, y = VolatileUbiquityBySample)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_line(
    aes(x = HarvestDate, y = fit.harv_dt_vs_num_vols$fitted.values),  # predicted data
    color = 'black') +
  geom_text(size=8,aes(label = sprintf("R² = %f", summary(fit.harv_dt_vs_num_vols)$r.squared), x = 280, y = 66)) +
  geom_text(size=8,aes(label = sprintf("p-value = %s", format(summary(fit.harv_dt_vs_num_vols)$coefficients[2, 4], scientific = TRUE)), x = 280, y = 64)) +
  GLOBAL_THEME +
  theme(
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24)
  ) +
  xlab("Harvest Date (julian days)") +
  ylab("Number of volatiles\n detected")
ggsave(
  "analyses/pca/figures/fig_3c.png",
  plot = fig_3c.plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px"
)

# Figure 3D
fit.harv_dt_vs_tva <- lm(pcs.df$TotalVolatileAbundanceBySample ~ pcs.df$HarvestDate)
summary(fit.harv_dt_vs_tva)
fig_3d.plot <- ggplot(pcs.df, aes(x = HarvestDate, y = TotalVolatileAbundanceBySample)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_line(
    aes(x = HarvestDate, y = fit.harv_dt_vs_tva$fitted.values),  # predicted data
    color = 'black') +
  geom_text(size=8, aes(label = sprintf("R² = %f", summary(fit.harv_dt_vs_tva)$r.squared), x = 280, y = 750)) +
  geom_text(size=8, aes(label = sprintf("p-value = %s", format(summary(fit.harv_dt_vs_tva)$coefficients[2, 4], scientific = TRUE)), x = 280, y = 720)) +
  GLOBAL_THEME +
  theme(
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24)
  ) +
  xlab("Harvest Date (julian days)") +
  ylab("Total volatile abundance (TIC)")
ggsave(
  "analyses/pca/figures/fig_3d.png",
  plot = fig_3d.plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px"
)

fig_3.plot <- ggarrange(
  fig_3a.plot,
  fig_3b.plot,
  fig_3c.plot,
  fig_3d.plot,
  nrow = 2,
  ncol = 2,
  labels = c("A", "B", "C", "D")
)
ggsave(
  "figures/Figure_3.png",
  plot = fig_3.plot,
  bg = "white",
  width = 10,
  height = 10,
  units = "in"
)

fig_4.df <- data.frame(
  ButylAcetate = pheno_tbl$`Butyl acetate`,
  HexylAcetate = pheno_tbl$`Hexyl acetate`,
  `2-methylbutyl acetate` = pheno_tbl$`2-methylbutyl acetate`,
  Hexanal = pheno_tbl$Hexanal,
  EthylButyrate = pheno_tbl$`Ethyl butyrate`,
  `1-Hexanol` = pheno_tbl$`1-Hexanol`,
  HarvestDate = pheno_tbl$date_jul_17_harv
)
formatted_names <- c(
  "Butyl acetate",
  "Hexyl acetate",
  "2-methylbutyl acetate",
  "Hexanal",
  "Ethyl butyrate",
  "1-Hexanol"
)

p_r2_val_y_pos <- list(
  c(300, 280),
  c(0, 0),
  c(0, 0),
  c(0, 0),
  c(0, 0),
  c(0, 0)
)

plots <- list()
for (i in 1:(length(fig_4.df) - 1)) {

  name <- names(fig_4.df)[i]
  formatted_name <- formatted_names[i]

  pheno_var <- fig_4.df[, name]
  hv_dt <- fig_4.df[, "HarvestDate"]

  y1 <- p_r2_val_y_pos[[i]][1]
  y2 <- p_r2_val_y_pos[[i]][2]

  print(name)
  print(formatted_name)
  print(paste0("Y1: ", y1, " -> Y2: ", y2))

  dat.df <- data.frame(Pheno = pheno_var, HarvestDate = hv_dt)

  # fit the model
  fit <- lm(dat.df$Pheno ~ dat.df$HarvestDate)
  plots[[i]] <- ggplot(dat.df, aes(x = HarvestDate, y = Pheno)) +
    geom_point(alpha = 0.4, size = 2) +
    geom_line(
      aes(x = HarvestDate, y = fit$fitted.values),  # predicted data
      color = 'black', alpha = 0.3) +
    geom_text(aes(label = sprintf("R² = %f", summary(fit)$r.squared), x = 280, y = y1)) +
    geom_text(aes(label = sprintf("p-value = %s", format(summary(fit)$coefficients[2, 4], scientific = TRUE)), x = 280, y = y2)) +
    GLOBAL_THEME +
    xlab("Harvest Date (julian days)") +
    ylab(paste0(formatted_name))
}
fig_4.plot <- ggarrange(plotlist = plots, ncol = 2, nrow = 3)
ggsave(
  filename = "analyses/pca/figures/fig_4.png",
  plot = fig_4.plot,
  bg = "white",
  width = 10,
  height = 10,
  units = "in"

)


# supplementary figure 1
# A
# total_sample_ubiquity_by_volatiles_plot.png
# B

######################################
## PCs TO ABC PHENOTYPE CORRELATION ##
######################################

abc_phenotypes <- pheno_tbl[, c(1, 112:ncol(pheno_tbl))]
dim(abc_phenotypes)
# [1] 515 42

abc_phenotypes.noaid <- abc_phenotypes[, -1]
abc_phenotypes.noaid <- as.matrix(abc_phenotypes.noaid)
dim(abc_phenotypes.noaid)
# [1] 515 41

# function for getting pcs
pcs <- function(x)
  pca_data.pca$x[, x]

# calculate the linear regression r2 value for PC1, PC2, PC3 and PC5 with ABC phenotypes
structure_matrix <- matrix(NA, nrow = ncol(abc_phenotypes.noaid), ncol = 5)
for (i in seq_len(ncol(abc_phenotypes.noaid))) {
  # pc1 linear model
  lm1_r2 <- summary(lm(abc_phenotypes.noaid[, i] ~ pcs(1)))$r.squared
  structure_matrix[i, 1] <- lm1_r2

  # pc1 and 2 linear model
  lm2_r2 <-
    summary(lm(abc_phenotypes.noaid[, i] ~ pcs(1) + pcs(2)))$r.squared - lm1_r2
  structure_matrix[i, 2] <- lm2_r2

  # pc1, pc2 and pc3 linear model
  lm3_r2 <-
    summary(lm(abc_phenotypes.noaid[, i] ~ pcs(1) + pcs(2) + pcs(3)))$r.squared - (lm2_r2 + lm1_r2)
  structure_matrix[i, 3] <- lm3_r2

  # pc1, pc2, pc3, pc4 linear model
  lm4_r2 <-
    summary(lm(abc_phenotypes.noaid[, i] ~ pcs(1) + pcs(2) + pcs(3) + pcs(4)))$r.squared - (lm3_r2 + lm2_r2 + lm1_r2)
  structure_matrix[i, 4] <- lm4_r2

  # pc1, pc2, pc3, pc4, pc5 linear model
  lm5_r2 <-
    summary(lm(abc_phenotypes.noaid[, i] ~ pcs(1) + pcs(2) + pcs(3) + pcs(4) + pcs(5)))$r.squared - (lm4_r2 + lm3_r2 + lm2_r2 + lm1_r2)
  structure_matrix[i, 5] <- lm5_r2

}
structure_matrix <- structure_matrix * 100


# get the minimum and maximum of variance
summary(rowSums(structure_matrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5659  3.1355  4.6828  9.3754 12.2012 31.8838

rownames(structure_matrix) <- colnames(abc_phenotypes)[-1]

# laod the table which contains the ordering of ABC phenotypes
phenotype_ordering_tbl <- read_excel('data/raw/abc_phenotypes_ordered.xlsx')

# order the structure_matrix based on phenotype grouping
structure_matrix <-
  structure_matrix[match(phenotype_ordering_tbl$phenotypes, rownames(structure_matrix)),]
rownames(structure_matrix) <- phenotype_ordering_tbl$cleaned

dat <- reshape2::melt(structure_matrix)
colnames(dat) <- c("Phenotype", "Legend", "Variance")
dat[which(dat$Legend == 1), 'Legend'] <- "PC1"
dat[which(dat$Legend == 2), 'Legend'] <- "PC2"
dat[which(dat$Legend == 3), 'Legend'] <- "PC3"
dat[which(dat$Legend == 4), 'Legend'] <- "PC4"
dat[which(dat$Legend == 5), 'Legend'] <- "PC5"

variance_plot <- generate_variance_plot(dat)
ggsave(
  filename = "analyses/pca/figures/pca_variance.png",
  plot = variance_plot,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in"
)

# Only keep selected phenotypes
dat_2017 <- dat[!grepl("2016", dat$Phenotype),]
dat_2017 <- dat_2017[-grep("after storage", dat_2017$Phenotype),]
dat_2017 <- dat_2017[-which(dat_2017$Phenotype == "FRAP"),]
dat_2017 <- dat_2017[-which(dat_2017$Phenotype == "Ripening time (2017)"),]
dat_2017 <- dat_2017[-which(dat_2017$Phenotype == "Δ Bx/acid (2017)"),]

# remove "at harvest" from the phenotype names"
dat_2017$Phenotype <- sub("at harvest ", "", dat_2017$Phenotype)

variance_plot_2017 <- generate_variance_plot(
  dat_2017,
  options =  list(
    axis_title_size = 25,
    axis_text_size = 25,
    legend_text_size = 20
  )
)
ggsave(
  filename = "analyses/pca/figures/pca_variance_subset.png",
  plot = variance_plot_2017,
  dpi = 600,
  width = 20,
  height = 8,
  units = "in"
)