################
## FIGS INDEX ##
################

## fig_3a - PCA bi-plot with Harvest Date
## fig_3b - Scatter plot of PC1 and harvest date with linear correlation
## fig_3c - Scatter plot of the ‘number of volatiles detected’ and harvest date of each sample with linear correlation
## fig_3d - Scatter plot of ‘total volatile abundance’ and harvest date of each sample with linear correlation

############################
## IMPORTS & DATA LOADING ##
############################

library(readxl)
library(tidyverse)
library(openxlsx)
library(ggpubr)
source("analyses/pca/utils.R")
source("themes/theme_main.R")
source("utils/basic_stats.R")
# library(viridis)


# loading the GCMS phenotype table
gcms_pheno_tbl <-
  read_excel("data/processed/Supplementary_Data.xlsx", sheet = "GCMS Data")
dim(gcms_pheno_tbl)
# [1] 515 107

# loading the ABC phenotype table
abc_pheno_tbl <- read_excel(
    "data/processed/Supplementary_Data.xlsx",
    sheet = "Apple Traits Data"
  )
dim(abc_pheno_tbl)
# [1] 515 44

# merge the GCMS and ABC phenotype data
pheno_tbl <- left_join(
  gcms_pheno_tbl, abc_pheno_tbl, by = c("appleid" = "apple_id")
)
dim(pheno_tbl)
# [1] 515 150

# get the column indices for aroma compounds (minus the appleid column)
aroma_columns_idx <- seq_along(gcms_pheno_tbl)[-1]

# Gather the PCA data
pca_data <- pheno_tbl[, aroma_columns_idx]
dim(pca_data)
# [1] 515 106

##################
## PCA ANALYSIS ##
##################

# generate the pca data frame
pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)
summary(pca_res)

# calculate the basic aroma stats like abundance and ubiquity
aroma_stats <- get_aroma_stats_by_samples(pca_data)


# calculate the proportion of variance of PCs
pov <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

pcs_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  AppleId = pheno_tbl$appleid,
  Use = pheno_tbl$use,
  World = pheno_tbl$world,
  HarvestDate = pheno_tbl$date_jul_17_harv,
  VolatileUbiquityBySample = aroma_stats$Ubiquity,
  TotalVolatileAbundanceBySample = aroma_stats$Abundance
)

# write the PCA data to a file
openxlsx::write.xlsx(
  pcs_df,
  file  = "analyses/pca/data/pcs_df.xlsx",
)

lim <- c(
  min(pcs_df$PC1), max(pcs_df$PC1),
  min(pcs_df$PC2), max(pcs_df$PC2)
)
limits <- c(min(lim), max(lim))

########################
## GENERATE PCA PLOTS ##
########################

# pca plot for harvest date change
fig_3a_plot <- generate_pca_biplot(
  pcs_df,
  choices = c("PC1", "PC2"),
  color_phenotype = c("HarvestDate", "Harvest Date\n (julian days)"),
  limits = limits,
  proportion_of_variance = pov
)
ggsave(
  "analyses/pca/figs/fig_3a.pdf",
  fig_3a_plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px",
  device = cairo_pdf
)

# PCA bi-plot for powerpoint presentation
fig_3a_plot_ppt <- generate_pca_biplot(
  pcs_df,
  choices=c("PC1", "PC2"),
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
  "analyses/pca/figs/fig_3a_ppt.png",
  fig_3a_plot_ppt,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px"
)

# fig_3b
fit_harv_dt_vs_pc1 <- lm(pcs_df$PC1 ~ pcs_df$HarvestDate)
summary(fit_harv_dt_vs_pc1)
fig_3b_plot <- ggplot(pcs_df, aes(x = HarvestDate, y = PC1)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_line(
    aes(x = HarvestDate, y = fit_harv_dt_vs_pc1$fitted.values),
    color = "black") +
  annotate(
    "text", size = 8, label = sprintf("R² = %0.2f", 
    summary(fit_harv_dt_vs_pc1)$r.squared), x = 280, y = -15
  ) +
  annotate("text", size = 8, label = sprintf("p = %0.2e",
    summary(fit_harv_dt_vs_pc1)$coefficients[2, 4]), x = 280, y = -17) +
  theme_classic2() +
  theme(
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  ) +
  xlab("Harvest Date (julian days)") +
  ylab("PC1")
ggsave(
  "analyses/pca/figs/fig_3b.pdf",
  plot = fig_3b_plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px",
  device = cairo_pdf
)

# Figure 3C - Scatter plot of the ‘number of volatiles detected’ and harvest date of each sample with linear correlation
fit_harv_dt_vs_num_vols <-
  lm(pcs_df$VolatileUbiquityBySample ~ pcs_df$HarvestDate)
summary(fit_harv_dt_vs_num_vols)
fig_3c_plot <- ggplot(
    pcs_df, aes(x = HarvestDate, y = VolatileUbiquityBySample)
  ) +
  geom_point(alpha = 0.4, size = 2) +
  geom_line(
    aes(x = HarvestDate, y = fit_harv_dt_vs_num_vols$fitted.values),
    color = "black") +
  annotate("text", size = 8, label = sprintf("R² = %0.2f",
    summary(fit_harv_dt_vs_num_vols)$r.squared), x = 280, y = 66) +
  annotate("text", size=8, label = sprintf("p-value = %0.2e",
    summary(fit_harv_dt_vs_num_vols)$coefficients[2, 4]), x = 280, y = 64) +
  theme_classic2() +
  theme(
    axis.title.x = element_text(size=11),
    axis.title.y = element_text(size=11)
  ) +
  xlab("Harvest Date (julian days)") +
  ylab("Number of volatiles detected")
ggsave(
  "analyses/pca/figs/fig_3c.pdf",
  plot = fig_3c_plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px",
  device = cairo_pdf
)

# Figure 3D - Scatter plot of ‘total volatile abundance’ and harvest date of each sample with linear correlation
fit_harv_dt_vs_tva <-
  lm(pcs_df$TotalVolatileAbundanceBySample ~ pcs_df$HarvestDate)
summary(fit_harv_dt_vs_tva)
fig_3d_plot <- ggplot(
    pcs_df, aes(x = HarvestDate, y = TotalVolatileAbundanceBySample)
  ) +
  geom_point(alpha = 0.4, size = 2) +
  geom_line(
    aes(x = HarvestDate, y = fit_harv_dt_vs_tva$fitted.values),
    color = "black") +
  annotate("text", size = 8, label = sprintf("R² = %0.2f",
    summary(fit_harv_dt_vs_tva)$r.squared), x=280, y=750
  ) +
  annotate("text", size = 8, label = sprintf("p-value = %0.2e",
    summary(fit_harv_dt_vs_tva)$coefficients[2, 4]), x=280, y=720
  ) +
  theme_classic2() +
  theme(
    axis.title.x = element_text(size=11),
    axis.title.y = element_text(size=11)
  ) +
  xlab("Harvest Date (julian days)") +
  ylab("Total volatile abundance (TIC)")
ggsave(
  "analyses/pca/figs/fig_3d.pdf",
  plot = fig_3d_plot,
  bg = "white",
  width = 815,
  height = 515,
  dpi = 75,
  units = "px",
  device = cairo_pdf
)

fig_3_plot <- ggarrange(
  fig_3a_plot,
  fig_3b_plot,
  fig_3c_plot,
  fig_3d_plot,
  nrow = 2,
  ncol = 2,
  labels = c("A", "B", "C", "D")
)
ggsave(
  "analyses/pca/figs/fig_3.pdf",
  plot = fig_3_plot,
  bg = "white",
  width = 15,
  height = 10,
  units = "in",
  device = cairo_pdf
)


######################################
## PCs TO ABC PHENOTYPE CORRELATION ##
######################################

abc_pheno_tbl <- pheno_tbl[, c(1, 112:ncol(pheno_tbl))]
dim(abc_pheno_tbl)
# [1] 515 40

abc_pheno_noaid_tbl <- abc_pheno_tbl[, -1]
abc_pheno_noaid_tbl <- as.matrix(abc_pheno_noaid_tbl)
dim(abc_pheno_noaid_tbl)
# [1] 515 39

# function for getting pcs
pcs <- function(x) {
  pca_res$x[, x]
}

# calculate the linear regression r2 value for PC1, PC2, PC3 and PC5 with ABC
#  phenotypes
structure_matrix <- matrix(NA, nrow = ncol(abc_pheno_noaid_tbl), ncol = 5)
for (i in seq_len(ncol(abc_pheno_noaid_tbl))) {
  # pc1 linear model
  lm1_r2 <- summary(lm(abc_pheno_noaid_tbl[, i] ~ pcs(1)))$r.squared
  structure_matrix[i, 1] <- lm1_r2

  # pc1 and 2 linear model
  lm2_r2 <-
    summary(lm(abc_pheno_noaid_tbl[, i] ~ pcs(1) + pcs(2)))$r.squared - lm1_r2
  structure_matrix[i, 2] <- lm2_r2

  # pc1, pc2 and pc3 linear model
  lm3_r2 <-
    summary(lm(abc_pheno_noaid_tbl[, i] ~ pcs(1) + pcs(2) + pcs(3)))$r.squared - (lm2_r2 + lm1_r2)
  structure_matrix[i, 3] <- lm3_r2

  # pc1, pc2, pc3, pc4 linear model
  lm4_r2 <-
    summary(lm(abc_pheno_noaid_tbl[, i] ~ pcs(1) + pcs(2) + pcs(3) + pcs(4)))$r.squared - (lm3_r2 + lm2_r2 + lm1_r2)
  structure_matrix[i, 4] <- lm4_r2

  # pc1, pc2, pc3, pc4, pc5 linear model
  lm5_r2 <-
    summary(lm(abc_pheno_noaid_tbl[, i] ~ pcs(1) + pcs(2) + pcs(3) + pcs(4) + pcs(5)))$r.squared - (lm4_r2 + lm3_r2 + lm2_r2 + lm1_r2)
  structure_matrix[i, 5] <- lm5_r2

}
structure_matrix <- structure_matrix * 100


# get the minimum and maximum of variance
summary(rowSums(structure_matrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5659  3.1355  4.6828  9.3754 12.2012 31.8838

rownames(structure_matrix) <- colnames(abc_pheno_tbl)[-1]

# laod the table which contains the ordering of ABC phenotypes
phenotype_ordering_tbl <- read_excel("data/raw/abc_phenotypes_ordered.xlsx")

# order the structure_matrix based on phenotype grouping
structure_matrix <-
  structure_matrix[match(phenotype_ordering_tbl$phenotypes, rownames(structure_matrix)),]
rownames(structure_matrix) <- phenotype_ordering_tbl$cleaned

dat <- reshape2::melt(structure_matrix)
colnames(dat) <- c("Phenotype", "Legend", "Variance")
dat[which(dat$Legend == 1), "Legend"] <- "PC1"
dat[which(dat$Legend == 2), "Legend"] <- "PC2"
dat[which(dat$Legend == 3), "Legend"] <- "PC3"
dat[which(dat$Legend == 4), "Legend"] <- "PC4"
dat[which(dat$Legend == 5), "Legend"] <- "PC5"

variance_plot <- generate_variance_plot(dat)
ggsave(
  filename = "analyses/pca/figs/pca_variance.pdf",
  plot = variance_plot,
  dpi = 600,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
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
  filename = "analyses/pca/figs/pca_variance_subset.pdf",
  plot = variance_plot_2017,
  dpi = 600,
  width = 20,
  height = 8,
  units = "in",
  device = cairo_pdf
)
