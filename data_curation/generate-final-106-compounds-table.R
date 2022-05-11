# Title     : Final curated compounds
# Objective : This script generates a final list of compounds that are used
#             for downstream analysis after they have been manually curated
#             based on the comments from our AAFC collaborators
# Created by: tayabsoomro
# Created on: 2021-06-09

library(readxl)

# load the final compounds data received from AAFC collaborators
fnl_cmpds <- read_excel('data/processed/final_compounds_table_rev2.xlsx')

cleaned.df <- read.table(
  'data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated_109.tsv',
  header = TRUE
)

# remove Pentane, 1-chloro-
cleaned.df <- cleaned.df[-which(cleaned.df$Name == 'Pentane, 1-chloro-'),]
fnl_cmpds <- fnl_cmpds[-which(fnl_cmpds$Name == 'Pentane, 1-chloro-'),]

# remove Methylene-chloride
cleaned.df <- cleaned.df[-grep("Methylene chloride", cleaned.df$Name),]
fnl_cmpds <- fnl_cmpds[-which(fnl_cmpds$Name == 'Methylene-chloride'),]

# remove 2-4-Hexadienal-294
cleaned.df <- cleaned.df[-grep("2,4-Hexadienal-294",cleaned.df$Name),]
fnl_cmpds <- fnl_cmpds[-which(fnl_cmpds$Name == '2-4-Hexadienal-294'),]

# renaming 2-4-Hexadienal-430 to be just 2-4-Hexadienal as it is the correct one
cleaned.df[grep("2,4-Hexadienal",cleaned.df$Name),'Name'] <- '2,4-Hexadienal'
fnl_cmpds[which(fnl_cmpds$Name == '2-4-Hexadienal-430'),'Name'] <- '2-4-Hexadienal'

# renaming 5-9-Undecadien-2-one, 6-10-dimethyl to (E)-Geranylacetone
cleaned.df[which(cleaned.df$Name == "5,9-Undecadien-2-one, 6,10-dimethyl"),'Name'] <- "(E)-Geranylacetone"
fnl_cmpds[which(fnl_cmpds$Name == "5-9-Undecadien-2-one, 6-10-dimethyl"),'Name'] <- "(E)-Geranylacetone"

# renaming 4-Oxohex-2-enal-571 to 4-Oxohex-2-enal
cleaned.df[which(cleaned.df$Name == "4-Oxohex-2-enal-571"),'Name'] <- "4-Oxohex-2-enal"
fnl_cmpds[which(fnl_cmpds$Name == '4-Oxohex-2-enal-571'),'Name'] <- '4-Oxohex-2-enal'

# renaming 4-Oxohex-2-enal-684 to 5-ethyl-2(5H)-furanone
cleaned.df[which(cleaned.df$Name == "4-Oxohex-2-enal-684"),'Name'] <- "5-ethyl-2(5H)-furanone"
fnl_cmpds[which(fnl_cmpds$Name == '4-Oxohex-2-enal-684'),'Name'] <- '5-ethyl-2(5H)-furanone'

# renaming 2-Buten-1-one, 1, 2-6-6-trimethyl-1-3-cyclohexadien-1-yl, to β-Damascenone
cleaned.df[which(cleaned.df$Name == "2-Buten-1-one, 1-(2,6,6-trimethyl-1,3-cyclohexadien-1-yl)-"),'Name'] <- "β-Damascenone"
fnl_cmpds[which(fnl_cmpds$Name == '2-Buten-1-one, 1, 2-6-6-trimethyl-1-3-cyclohexadien-1-yl,'),'Name'] <- 'β-Damascenone'

# renaming 3-Hexen-1-ol, acetate to 3-Hexenyl acetate
cleaned.df[which(cleaned.df$Name == "3-Hexen-1-ol, acetate"),'Name'] <- "3-Hexenyl acetate"
fnl_cmpds[which(fnl_cmpds$Name == '3-Hexen-1-ol, acetate'),'Name'] <- '3-Hexenyl acetate'

# renaming 3-Methyl-3-buten-1-ol, acetate to 3-Methyl-3-butenyl acetate
cleaned.df[which(cleaned.df$Name == "3-Methyl-3-buten-1-ol, acetate"),'Name'] <- "3-Methyl-3-butenyl acetate"
fnl_cmpds[which(fnl_cmpds$Name == '3-Methyl-3-buten-1-ol, acetate'),'Name'] <- '3-Methyl-3-butenyl acetate'

# renaming 5-Hexene-1-ol, acetate to 5-Hexenyl acetate
cleaned.df[which(cleaned.df$Name == "5-Hexene-1-ol, acetate"),'Name'] <- "5-Hexenyl acetate"
fnl_cmpds[which(fnl_cmpds$Name == '5-Hexene-1-ol, acetate'),'Name'] <- '5-Hexenyl acetate'

#######################################
## MISSINGNESS THRESHOLD FILTERATION ##
#######################################

# seeing the change in number of compounds for various missingness thresholds
pmissing <- seq(0,100,5)
num_compounds <- NULL
for(p in pmissing){
  num_compounds <- c(num_compounds,nrow(fnl_cmpds[fnl_cmpds$`Samples (%)` > (100-p),]))
}
plot(
  num_compounds,
  pmissing,
  main = "Number of compounds retained for varying sample \nmissingness percent", 
  xlab = "Number of compounds", ylab = "Sample Missingness (%)"
)

# Generating a final cleaned GCMS data table
write.table(
  fnl_cmpds,
  'data/processed/final_compounds_table_rev2_addressed.tsv',
  sep = "\t",
  row.names = FALSE
)

# write the cleaned GCMS table
write.table(
  cleaned.df,
  'data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated_106.tsv',
  sep = "\t",
  row.names = FALSE
)

