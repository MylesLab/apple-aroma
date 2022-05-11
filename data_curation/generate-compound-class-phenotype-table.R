library(readxl)
library(xlsx)
library(tidyverse)


source('data_curation/generate-compound-class-corrected-gcms-table/util.R')

# load the GCMS phenotype table
gcms_pheno_tbl <- read_excel('data/processed/Supplementary_Data.xlsx', sheet = 'GCMS Data')
dim(gcms_pheno_tbl)
#[1] 515 107

# get all the compound classes
all_classes <- unique(classification_pivot$Classification)

# create a matrix to store the new class GCMS data
class_gcms_pheno_tbl <- matrix(, nrow = nrow(gcms_pheno_tbl), ncol = length(all_classes))
dim(class_gcms_pheno_tbl)
# [1] 515 13

for (sample in seq_len(nrow(gcms_pheno_tbl))) {
  colIdx <- 1
  for (class in all_classes) {
    compounds <- get_all_compounds_with_class(class)

    class_gcms_pheno_tbl[sample, colIdx] <- rowSums(gcms_pheno_tbl[sample, compounds])

    colIdx <- colIdx + 1
  }
}

class_gcms_pheno_tbl.df <- as.data.frame(class_gcms_pheno_tbl)
colnames(class_gcms_pheno_tbl.df) <- all_classes

# add appleid column at the start of the class_gcms_pheno_tbl.df
class_gcms_pheno_tbl.df <- class_gcms_pheno_tbl.df %>%
  add_column(appleid = gcms_pheno_tbl$appleid, .before = 1)


# write the each class data into individual file for GWAS
for (column in colnames(class_gcms_pheno_tbl.df)[-1]) {
  file_name <- paste0(
    gsub("/", "-", gsub(" ", "_", column)),
    ".txt"
  )
  print(file_name)
  write.table(
    class_gcms_pheno_tbl.df[, c('appleid', column)],
    file = paste0('data/processed/phenotype_tables/classes/', file_name),
    row.names = FALSE,
    col.names = F,
    sep = " "
  )
}

# write the class GCMS table to a file
write.xlsx(
  class_gcms_pheno_tbl.df,
  file = 'data/processed/Supplementary_Data.xlsx',
  row.names = FALSE,
  sheetName = 'Class GCMS Data',
  append = TRUE
)
