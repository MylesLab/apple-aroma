library(readxl)

gcms_class_corrected_pheno_tbl <- read_excel(
  'data/processed/Supplementary_Data.xlsx',
  sheet = "Class Corrected GCMS Data"
)
dim(gcms_class_corrected_pheno_tbl)
# [1] 515 107


# write the each compound data into individual file for GWAS
for (column in colnames(gcms_class_corrected_pheno_tbl)[-1]) {
  file_name <- paste0(
    gsub("/", "-", gsub(" ", "_", column)),
    ".txt"
  )
  write.table(
    gcms_class_corrected_pheno_tbl[, c('appleid', column)],
    file = paste0('data/processed/phenotype_tables/class_corrected/', file_name),
    row.names = FALSE,
    col.names = F,
    sep = " "
  )
}




