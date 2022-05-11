library(readxl)
library(dplyr)


gcms_pheno_tbl <- read_excel('data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx')
dim(gcms_pheno_tbl)
# [1] 515 107

# the order of apple ids that every other df should follow
order <- gcms_pheno_tbl$appleid

classification_pivot_tbl <- read_excel('data/processed/classification_pivot.xlsx')
dim(classification_pivot_tbl)
# [1] 106 2

all_alcohols <- unlist(classification_pivot_tbl[classification_pivot_tbl$Classification == "Alcohol",]$Name)
all_esters <- unlist(classification_pivot_tbl[grep("Ester",classification_pivot_tbl$Classification),]$Name)

abc_pheno_tbl <- read_excel('data/processed/sup_tbl_2-abc_phenotype_table_v2.xlsx')
dim(abc_pheno_tbl)
# [1] 515 44

# re-organize the phenotype table so that it follows the same order as the gcms_class_tbl
abc_pheno_tbl <- abc_pheno_tbl[match(order, abc_pheno_tbl$apple_id),]

correlation_df <- data.frame(
  apple_id = order,
  HarvestDate = abc_pheno_tbl$date_jul_17_harv,
  Alcohol = gcms_class_tbl$Alcohol,
  EsterBranched = gcms_class_tbl$`Ester, Branched chain`,
  EsterStraight = gcms_class_tbl$`Ester, Straight chain`
)