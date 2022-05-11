library(readxl)
library(tictoc)
library(hash)
library(dplyr)
library(tibble)
library(beepr)
library(xlsx)

source('data_curation/generate-compound-class-corrected-gcms-table/util.R')

# read the GCMS phenotype table
gcms_pheno_tbl <- read_excel('data/processed/Supplementary_Data.xlsx', sheet = 'GCMS Data')
dim(gcms_pheno_tbl)
# [1] 515 107

# GCMS table without apple id
gcms_pheno_tbl.noaid <- gcms_pheno_tbl[, 2:ncol(gcms_pheno_tbl)]
dim(gcms_pheno_tbl.noaid)
# [1] 515 106

# create a matrix to store the new compound class corrected data
new_gcms_pheno_tbl <- matrix(, nrow = nrow(gcms_pheno_tbl), ncol = ncol(gcms_pheno_tbl) - 1)
dim(new_gcms_pheno_tbl)

tictoc::tic("OUTER FOR LOOP")
for (sample in seq_len(nrow(gcms_pheno_tbl))) {
    class_total_abund_hash <- hash()
    for (compound in seq_along(gcms_pheno_tbl)[-1]) {
        compound_name <- colnames(gcms_pheno_tbl)[compound]
        compound_class <- get_compound_class(compound_name)
        compound_abundance <- as.numeric(gcms_pheno_tbl[sample, compound])

        # if the total abundance ofr this sample for this class is already calculated, then use that
        # other wise, calculate it
        if (has.key(compound_class, class_total_abund_hash)) {
            total_class_abundance <- class_total_abund_hash[[compound_class]]
        } else {
            compound_names_with_class <- get_all_compounds_with_class(compound_class)
            total_class_abundance <- rowSums(gcms_pheno_tbl[sample, compound_names_with_class])
            class_total_abund_hash[[compound_class]] <- total_class_abundance
        }

        new_val <- sum(
          ifelse(
            total_class_abundance == 0,
            0,
            compound_abundance / total_class_abundance
          )
        )

        new_gcms_pheno_tbl[sample, compound - 1] <- new_val
    }
}
tictoc::toc()
beep(sound = 3)

class_corrected_gcms_pheno_tbl <- as.data.frame(new_gcms_pheno_tbl)
colnames(class_corrected_gcms_pheno_tbl) <- colnames(gcms_pheno_tbl)[-1]

class_corrected_gcms_pheno_tbl <- class_corrected_gcms_pheno_tbl %>%
  add_column(appleid = gcms_pheno_tbl$appleid, .before = 1)

class_corrected_gcms_pheno_tbl.noaid <- class_corrected_gcms_pheno_tbl[, 2:ncol(class_corrected_gcms_pheno_tbl)]

dim(class_corrected_gcms_pheno_tbl)

sum(rowSums(class_corrected_gcms_pheno_tbl.noaid == 0))
sum(rowSums(gcms_pheno_tbl.noaid == 0))

# export the dataframe
write.xlsx(
  class_corrected_gcms_pheno_tbl,
  file = "data/processed/Supplementary_Data.xlsx",
  sheetName = "Class Corrected GCMS Data 2",
  row.names = FALSE,
  append = TRUE
)

test.df <- data.frame(
  Corrected = as.numeric(class_corrected_gcms_pheno_tbl.noaid[1,]),
  Uncorrected = as.numeric(gcms_pheno_tbl.noaid[1,])
)

ggplot(test.df, aes(Corrected, Uncorrected)) +
  stat_summary(fun.data = mean_cl_normal, alpha = 0.5) +
  geom_smooth(method = 'lm', show.legend = T, formula = y ~ x) +
  GLOBAL_THEME +
  xlab("Corrected") +
  ylab("Uncorrected")


