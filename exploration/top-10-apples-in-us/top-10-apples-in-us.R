library(xlsx)
library(ggplot2)

# loading the GCMS data
gcms_pheno_tbl <- read_excel('data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx')

# get the list of apple ids as base for sorting the metadata table
order <- gcms_pheno_tbl$appleid

# loading the apple ABC metadata
apple_metadata_tbl <- read.delim('data/raw/pheno_meta_data.csv', sep = ",")

# sort the rows so that they are similar to GCMS table
apple_metadata_tbl <- apple_metadata_tbl[match(order, apple_metadata_tbl$apple_id),]


# here are the top 10 US apple varieties
top_10_us_varieties_names <- c(
  "Red Delicious",
  "Gala",
  "Granny Smith",
  "Fuji Red Sport Type 2",
  "Golden Delicious",
  "Honeycrisp",
  "McIntosh",
  "Rome Beauty Law",
  "Pink Lady",
  "Empire"
)


for (var in top_10_us_varieties_names) {
  apple_id <- apple_metadata_tbl[grep(paste0("^", var, "$"), apple_metadata_tbl$PLANTID), 'apple_id']
  print(apple_id)
}

#############
## RUN PCA ##
#############

# create PCA data
pca_data <- gcms_pheno_tbl[, 2:ncol(gcms_pheno_tbl)]
dim(pca_data)
# [1] 515 106

# run PCA
pca_data.pca <- prcomp(pca_data, center = TRUE, scale. = TRUE)

pcs.df <- data.frame(
  PC1 = pca_data.pca$x[, 1],
  PC2 = pca_data.pca$x[, 2],
  AppleID = gcms_pheno_tbl$appleid
)

ggplot(pcs.df, aes(x = PC1, y = PC2)) + geom_point()




