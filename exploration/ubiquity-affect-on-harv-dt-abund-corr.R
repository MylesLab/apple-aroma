# TODO: Do this exercise again and generate your own insights
# TODO: Figure out 34% of positively correlated vs. 42% negatively correlated with harvest date

library(dplyr)
library(readxl)

gcms_pheno_tbl <- read_excel('data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx')
dim(gcms_pheno_tbl)
# [1] 515 107

aroma_columns_idx <- 2:ncol(gcms_pheno_tbl)

abc_pheno_tbl <- read_excel('data/processed/sup_tbl_2-abc_phenotype_table_v2.xlsx')
abc_pheno_tbl <- abc_pheno_tbl[, c("apple_id", "date_jul_17_harv")]
dim(abc_pheno_tbl)
#[1] 515 2

pheno_tbl <- inner_join(gcms_pheno_tbl, abc_pheno_tbl, by = c("appleid" = "apple_id"))
dim(pheno_tbl)
# [1] 515 108

num_cmpds <- NULL
r_vals <- NULL
p_vals <- NULL
for (i in aroma_columns_idx) {
  cols_idx <- i:length(aroma_columns_idx)
  abund <- rowSums(pheno_tbl[, cols_idx])
  harv_dt <- pheno_tbl$date_jul_17_harv

  num_cmpds[i - 1] <- length(cols_idx)
  r_vals[i - 1] <- cor.test(abund, harv_dt)$estimate
  p_vals[i - 1] <- cor.test(abund, harv_dt)$p.value
}

plot(num_cmpds, r_vals, type = "o")
plot(num_cmpds, p_vals, type = "o")


