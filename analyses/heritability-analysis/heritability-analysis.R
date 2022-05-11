# Title     : SNP Heritability Analysis
# Objective : This script performs the SNP heritability analysis
# Created on: 2021-07-27


#####################
## LIBARAY IMPORTS ##
#####################

library(readr)
library(readxl)
library(ggplot2)
library(RColorBrewer)
source('analyses/heritability-analysis/utils.R')
source('themes/theme_avenir.R')

#######################################################
## GENERATING THE PHENOTYPE HERITABILITY INPUT TABLE ##
#######################################################

# get the individual ids (apple ids) from the grm table
indv <- read_table2(
  'data/raw/snp-heritability/20210805_gcms_grm.grm.id',
  col_names = FALSE
)

# order of the apple ids
order <- indv$X1

# Read the GRM
apple_grm <- ReadGRMBin(
  'data/raw/snp-heritability/20210805_gcms_grm',
  AllN = F, size = 4
)

# GRM dataset
off_dat <- apple_grm[["off"]]
diag <- apple_grm[["diag"]]
hist(diag, breaks = 100)

# get the gcms data
gcms_data <- read_excel(
  'data/processed/Supplementary_Data.xlsx',
  sheet = "GCMS Data"
)
dim(gcms_data)
# [1] 515 107

# only apple ids that are in GCMS data
order_in <- order[which(order %in% gcms_data$appleid)]

# these apple ids (accessions) do not have the genetic data and thus should be
# removed from GCMS data for SNP heritability
missing_aids <- as.numeric(
  unlist(gcms_data[which(!(gcms_data$appleid %in% order_in)),'appleid'])
)

# only the accessions for which we have the genetic data
gcms_data_cleaned <- gcms_data[which(!(gcms_data$appleid %in% missing_aids)),]


# perform PCA
gcms_data_ordered <- gcms_data_cleaned[match(order_in, gcms_data_cleaned$appleid),]
pca_data_gcms <- gcms_data_ordered[,2:ncol(gcms_data_ordered)]
gcms.pca <- prcomp(pca_data_gcms, scale. = TRUE, center = TRUE)

# get positions of samples along PCs
pos_pcs <- as.data.frame(gcms.pca$x)[1:5]

# add family id column
gcms_data_ordered <- cbind(appleid=gcms_data_ordered$appleid, gcms_data_ordered)

pheno_heritability <- cbind(gcms_data_ordered, pos_pcs)
write.table(
  pheno_heritability,
  file = "data/raw/snp-heritability/pheno_heritability.phen",
  sep = "\t", row.names = FALSE, col.names = FALSE
)

######################################
## VISUALIZING THE SNP HERITABILITY ##
######################################

# load the SNP heritability data
snp_herit_tbl <- read.table(
  'data/raw/snp-heritability/gcms_snp_heritability.tsv',
  sep = "\t",
  header = TRUE
)

snp_herit_tbl_ordered <- snp_herit_tbl[order(snp_herit_tbl$Heritability, decreasing = TRUE),]

herit_plot <- ggplot(snp_herit_tbl, aes(x = Heritability, y = reorder(Phenotype, Heritability), fill = Heritability)) +
  geom_bar(stat = "identity", color = "black", width = 1, position = position_dodge()) +
  theme_classic() +
  theme_avenir() +
  theme(
    axis.text.y = element_text(size = 5, hjust = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5, size = 12),
    axis.title.y = element_text(hjust = 0.5, size = 12)
  ) +
  xlab("Heritability (r²)") +
  ylab("") +
  scale_fill_continuous() +
  geom_errorbar(aes(xmin = Heritability - StdErr, xmax = Heritability + StdErr), width = 0.5,
                position = position_dodge(.5))

ggsave(
  'figures/heritability/heritability_plot.png',
  herit_plot,
  dpi = 600
)

# remove PCs from the data
snp_herit_tbl.no_pcs <- snp_herit_tbl[-grep("PC", snp_herit_tbl$Phenotype),]

# adding missingness to heritability calculation
missingness <- NULL
idx <- 1
max_val <- 0
for (name in snp_herit_tbl.no_pcs$Phenotype) {
  clean_name <- trimws(name)
  if (clean_name %in% colnames(gcms_data)) {
    val <- as.numeric(colSums(gcms_data[, clean_name] == 0))
    missingness[idx] <- val

    if (val > max_val) {
      max_val <- val
    }
    idx <- idx + 1
  } else {
    print(clean_name)
  }

}
snp_herit_tbl.no_pcs$Missingness <- missingness

View(snp_herit_tbl.no_pcs)


ggplot(snp_herit_tbl.no_pcs, aes(x = Heritability, y = Missingness)) + geom_count(alpha = 0.5)

plot(
  snp_herit_tbl$Heritability,
  snp_herit_tbl$Missingness,
  xlab = "Heritability (r²)",
  ylab = "Samples missing",
  title(main = "Missingness and heritability")
)
