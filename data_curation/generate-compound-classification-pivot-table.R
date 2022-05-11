library(readxl)
library(xlsx)

# read the table
classification_tbl <- read.csv(
  'data/processed/final_compounds_table_rev2_addressed.tsv',
  sep = "\t", header = TRUE
)
classification_tbl <- classification_tbl[, c("Name", "Classification")]

# export the pivot table
openxlsx::write.xlsx(
  classification_tbl,
  file = "data/processed/classification_pivot.xlsx",
  row.names = FALSE
)