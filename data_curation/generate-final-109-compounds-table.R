# Title     : Generate Compound Summary Table
# Objective : This script generates a compound summary table which will be sent
#             to our AAFC collaborators to find out which ones are the most likely
#             to be found in apple.
# Created by: tayabsoomro
# Created on: 2021-05-26

pheno_tbl <- read.table(
  'data/processed/standardized/gcms_phenotype_table_std_w_benzaldehyde_sim_gt_850_109_cmpds.tsv',
  header = TRUE
)

# moving nursery ids from colum to rownames
nurseryids <- pheno_tbl$nurseryid
pheno_tbl <- pheno_tbl[, !names(pheno_tbl) %in% "nurseryid"]
rownames(pheno_tbl) <- nurseryids

# transposed phenotype table
t_pheno_tbl <- t(pheno_tbl)

final.df <- matrix(, nrow=0, ncol=4)
colnames(final.df) <- c("Name", "Samples (%)", "Classification", "Likely")
for (name in row.names(t_pheno_tbl)) {
  
  # the abundance of samples for this compound
  abund <- (sum(t_pheno_tbl[name,] != 0) / ncol(t_pheno_tbl))*100
  
  final.df <- rbind(
    final.df,
    c(name, abund, "","")
  )
}

final.df <- as.data.frame(final.df)
final.df$`Samples (%)` <- as.numeric(final.df$`Samples (%)`)

final.df$Name <- gsub("^X","",final.df$Name)

write.table(
  final.df,
  'data/processed/final_compound_table.tsv',
  sep = "\t",
  row.names =  FALSE
)


