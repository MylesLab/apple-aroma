# Title     : Final Supplementary Tables
# Objective : This script generates all the supplementary tables that are
#             required for the publication
# Created by: tayabsoomro
# Created on: 2021-06-22

############################
## IMPORTS & DATA LOADING ##
############################

library(xlsx)
library(dplyr)
library(tibble)
library(readxl)

# load functions
source("data_curation/generate-supplementary-tables/utils.R")

# load the final data for compounds and samples
pheno_tbl <- read.table(
  "data/processed/standardized/gcms_phenotype_table_std_w_benzaldehyde_sim_gt_850_106_cmpds.tsv",
  header = TRUE
)
dim(pheno_tbl)
# [1] 553 107

# get the classification data
class_data <- read_excel("data/processed/final_compounds_table_rev2_addressed.xlsx",)

# get the apple id and nursery id pivot data
pivot_tbl <- read.table("data/raw/nursery-id_apple-id_pivot.tsv")

abc_pop_info_tbl <- read_excel("data/raw/20200204_abc_pop_info.xlsx", col_types = "skip")

# get the abc phenotype data
metadata <- read.csv("data/raw/pheno_meta_data.csv", header = TRUE)


pheno_tbl <- as.data.frame(pheno_tbl)

##################################################
## SUPPLEMENATARY TABLE 1: GCMS PHENOTYPE TABLE ##
##################################################

# get the apple ids for the nursery ids
aids <- NULL
for (i in seq_along(pheno_tbl$nurseryid)) {
  aids[i] <-
    pivot_tbl[which(pivot_tbl$nursery_id == pheno_tbl[i, "nurseryid"]), "apple_id"]
}

# make the column apple id and remove the apple id row
pheno_tbl <- pheno_tbl[, 2:ncol(pheno_tbl)]
pheno_tbl <- cbind(appleid = aids, pheno_tbl)

# cleanup the compound names
colnames(pheno_tbl) <- cleanup_compound_names(pheno_tbl)

# saving the final GCMS phenotype table
write.xlsx(
  as.data.frame(pheno_tbl),
  "data/processed/sup_tbl_1-final_gcms_phenotype_table.xlsx",
  row.names = FALSE
)

#################################################
## SUPPLEMENATARY TABLE 2: ABC PHENOTYPE TABLE ##
#################################################

# filter for only relevant pheno types
abc_pheno_data <-
  metadata[, !names(metadata) %in% c("PLANTID", "species", "release_year", "ACCID")]

# only keep the accessions which are used in GCMS
aids_to_keep <- pheno_tbl$appleid
aids_to_keep_idx <- which(abc_pheno_data$apple_id %in% aids_to_keep)

nrow(abc_pheno_data)
# [1] 1119

abc_pheno_data <- abc_pheno_data[aids_to_keep_idx,]

nrow(abc_pheno_data)
# [1] 553

# saving the final phenotype table
write.xlsx(
  as.data.frame(abc_pheno_data),
  "data/processed/sup_tbl_2-abc_phenotype_table.xlsx",
  row.names = FALSE
)

#################################################
## LIST OF APPLES TO SEND FOR SENSORY ANALYSIS ##
#################################################

# add PI IDs
pi_dat <-
  unique(abc_pop_info_tbl[which(abc_pop_info_tbl$ACP == "PI"), c("ACNO", "apple_id")])
colnames(pi_dat) <- c("PI", "appleid")
pheno_tbl <- right_join(pi_dat, pheno_tbl, by = "appleid")


zero_items <- rep(0, nrow(pheno_tbl))

pheno_tbl <-
  add_column(pheno_tbl, tot_Alcohol = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl, tot_Lactone = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl, tot_Ester = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl, tot_Hydrocarbon = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl,
             tot_C13Norisoprenoid = zero_items,
             .after = "appleid"
  )
pheno_tbl <-
  add_column(pheno_tbl, tot_Sesquiterpene = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl, tot_Furan = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl, tot_Ketone = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl, tot_Aldehyde = zero_items, .after = "appleid")
pheno_tbl <-
  add_column(pheno_tbl, tot_Acid = zero_items, .after = "appleid")

# add cumulative sum of various classes
rows <- seq_len(nrow(pheno_tbl))
cols <- 3:ncol(pheno_tbl)

for (r in rows) {
  all_cls <- NULL
  for (c in cols) {
    abundance <- as.numeric(pheno_tbl[r, c])
    if (abundance != 0) {
      classes_present <-
        as.character(class_data[which(class_data$Name == colnames(pheno_tbl)[c]), "Classification"])
      all_cls <- c(all_cls, classes_present)
    }
  }

  print(paste0("R: ", r))

  # figure out the total classifications
  tot_Alcohol <- sum(all_cls == "Alcohol")
  tot_Lactone <- sum(all_cls == "Lactone")
  tot_Ester <- sum(grep("Ester", all_cls))
  tot_Hydrocarbon <- sum(all_cls == "Hydrocarbon")
  tot_C13Norisoprenoid <- sum(all_cls == "C13-norisoprenoid")
  tot_Sesquiterpene <- sum(all_cls == "Sesquiterpene")
  tot_Furan <- sum(all_cls == "Furan")
  tot_Ketone <- sum(all_cls == "Ketone")
  tot_Aldehyde <- sum(all_cls == "Aldehyde")
  tot_Acid <- sum(all_cls == "Acid")

  # assign the total classifications in the rows
  pheno_tbl[r, "tot_Alcohol"] <- as.numeric(tot_Alcohol)
  pheno_tbl[r, "tot_Lactone"] <- tot_Lactone
  pheno_tbl[r, "tot_Ester"] <- tot_Ester
  pheno_tbl[r, "tot_Hydrocarbon"] <- tot_Hydrocarbon
  pheno_tbl[r, "tot_C13Norisoprenoid"] <- tot_C13Norisoprenoid
  pheno_tbl[r, "tot_Sesquiterpene"] <- tot_Sesquiterpene
  pheno_tbl[r, "tot_Furan"] <- tot_Furan
  pheno_tbl[r, "tot_Ketone"] <- tot_Ketone
  pheno_tbl[r, "tot_Aldehyde"] <- tot_Aldehyde
  pheno_tbl[r, "tot_Acid"] <- tot_Acid
}

# esters
top20esters <-
  head(pheno_tbl[order(pheno_tbl$tot_Ester, decreasing = TRUE), c("PI", "appleid")], n = 20)
top20esters$comments <- "Top 20 Esters"
top20esters$rank <- 1
bot20esters <-
  tail(pheno_tbl[order(pheno_tbl$tot_Ester, decreasing = TRUE), c("PI", "appleid")], n = 20)
bot20esters$comments <- "Bottom 20 Esters"
bot20esters$rank <- 1

# aldehydes
top10aldehydes <-
  head(pheno_tbl[order(pheno_tbl$tot_Aldehyde, decreasing = TRUE), c("PI", "appleid")], n = 10)
top10aldehydes$comments <- "Top 10 Aldehydes"
top10aldehydes$rank <- 2
bot10aldehydes <-
  tail(pheno_tbl[order(pheno_tbl$tot_Aldehyde, decreasing = TRUE), c("PI", "appleid")], n = 10)
bot10aldehydes$comments <- "Bottom 10 Aldehydes"
bot10aldehydes$rank <- 2

# alcohols
top5alcohols <-
  head(pheno_tbl[order(pheno_tbl$tot_Alcohol, decreasing = TRUE), c("PI", "appleid")],
       n =
         5
  )
top5alcohols$comments <- "Top 5 Alcohols"
top5alcohols$rank <- 4
bot5alcohols <-
  tail(pheno_tbl[order(pheno_tbl$tot_Alcohol, decreasing = TRUE), c("PI", "appleid")],
       n =
         5
  )
bot5alcohols$comments <- "Bottom 5 Alchols"
bot5alcohols$rank <- 4

# ketones
top10ketones <-
  head(pheno_tbl[order(pheno_tbl$tot_Ketone, decreasing = TRUE), c("PI", "appleid")],
       n =
         10
  )
top10ketones$comments <- "Top 10 Ketones"
top10ketones$rank <- 3
bot10ketones <-
  head(pheno_tbl[order(pheno_tbl$tot_Ketone, decreasing = TRUE), c("PI", "appleid")],
       n =
         10
  )
bot10ketones$comments <- "Bottom 10 Ketones"
bot10ketones$rank <- 3

# combine the sensory data together
sensory_data <- rbind(
  top20esters,
  bot20esters,
  top10aldehydes,
  bot10aldehydes,
  top10ketones,
  bot10ketones,
  top5alcohols,
  bot5alcohols
)

# get only the unique apple ids
sensory_data <-
  as.data.frame(sensory_data[!duplicated(sensory_data$appleid),])

# update the column names to match the appleid column between the two tables
colnames(sensory_data) <- c("PI", "apple_id", "comments", "rank")

# connect nursery data with sensory data
nursery_data <- read.table("data/raw/20161024_ABC_design_grid_zm.txt",
                           header = TRUE
)

# only keep the relevant columns in nursery id table
nursery_data <-
  nursery_data[, c("Row", "Tree", "Block", "nursery_id", "PLANTID", "apple_id")]

# join the nursery data with the sensory panel data
final_sensory_data_tbl <-
  right_join(nursery_data, sensory_data, by = "apple_id")

nrow(final_sensory_data_tbl)
# [1] 149

# get the phenotype metadata for attaching the harvest dates
metadata <- read.csv("data/raw/pheno_meta_data.csv",
                     header = TRUE
)

# obtain only the harvest dates with the apple id
dates <-
  metadata[, c("apple_id", "date_jul_16_harv", "date_jul_17_harv")]

# join the harvest dates to the sensory data
final_sensory_data_tbl <- left_join(final_sensory_data_tbl, dates)

# remove PI column from for the sensory data
final_sensory_data_tbl <-
  final_sensory_data_tbl[, !names(final_sensory_data_tbl) %in% "PI"]

# write the table to file
write.table(
  final_sensory_data_tbl,
  "data/processed/apples_shortlisted_for_sensory_panel_analysis.tsv",
  sep = "\t",
  row.names = FALSE
)

################################################
## LIST OF TOP ESTERS, ALDEHYDES AND ALCOHOLS ##
################################################

all_esters <-
  class_data[grep("Ester", class_data$Classification), c("Name", "Samples (%)")]
topesters <-
  as.data.frame(head(all_esters[order(all_esters$`Samples (%)`, decreasing = TRUE),], n = 20)[, c("Name")])
write.xlsx(
  topesters,
  file = "data/processed/classification_of_compounds.xlsx",
  sheetName = "Top20Esters",
  append = FALSE,
  row.names = FALSE
)

all_aldehydes <-
  class_data[grep("Aldehyde", class_data$Classification), c("Name", "Samples (%)")]
topaldehydes <-
  as.data.frame(head(all_aldehydes[order(all_aldehydes$`Samples (%)`, decreasing = TRUE),], n = 10)[, c("Name")])
write.xlsx(
  topaldehydes,
  file = "data/processed/classification_of_compounds.xlsx",
  sheetName = "Top10Aldehydes",
  append = TRUE,
  row.names = FALSE
)

all_alcohols <-
  class_data[grep("Alcohol", class_data$Classification), c("Name", "Samples (%)", "Classification")]
topalcohols <-
  as.data.frame(head(all_alcohols[order(all_alcohols$`Samples (%)`, decreasing = TRUE),], n = 5)[, c("Name")])
write.xlsx(
  topalcohols,
  file = "data/processed/classification_of_compounds.xlsx",
  sheetName = "Top5Alcohols",
  append = TRUE,
  row.names = FALSE
)

all_ketones <-
  class_data[grep("Ketone", class_data$Classification), c("Name", "Samples (%)")]
topketones <-
  as.data.frame(head(all_ketones[order(all_ketones$`Samples (%)`, decreasing = TRUE),], n = 5)[, c("Name")])
write.xlsx(
  topketones,
  file = "data/processed/classification_of_compounds.xlsx",
  sheetName = "Top5Ketones",
  append = TRUE,
  row.names = FALSE
)