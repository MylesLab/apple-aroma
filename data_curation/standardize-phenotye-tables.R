# Title     : Standardizing the phenotypic table
# Objective : This script generates the 3 standardized tables for each of the
#             standards
# Created by: tayabsoomro
# Created on: 2021-03-22

#######################
## TABLE OF CONTENTS ##
#######################

## 0. IMPORTS & DATA LOADING
## 1. DISTRIBUTION OF SAMPLE COVERAGE FOR INT. STDs.
## 2. GENERATING FINAL STANDARDIZATION DATAFRAME
## 3. GENERIC STANDARDIZATION ROUTINE
## 4. STANDARDIZE WITH 2-hexanone-1,1,1,3.3-d5
## 5. STANDARDIZE WITH Benzaldehyde-d6
## 6. STANDARDIZE WITH Ethyl acetate-d8

## 0. IMPORTS & DATA LOADING
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# importing the standards data frame
std.df <- read_excel(
  'data/raw/InternalStandardData.xlsx',
  col_types = "text"
)

#######################################################
## 1. DISTRIBUTION OF SAMPLE COVERAGE FOR INT. STDs. ##
#######################################################

# we need to figure out the distribution of number of samples for each of the
# internal standards

# get indexes for all the rows which have non-empty Area values
non_na_std_idxs <- which(!is.na(std.df$Area))

# create a table with name of standards, their area values and the nursery
# ids column
stds_w_area_values <- std.df[non_na_std_idxs, c('Name', 'Area', 'eUnit')]

nrow(stds_w_area_values)
# [1] 1792

# remove all the Blank rows
stds_w_area_values <- stds_w_area_values[-grep(".*Blank.*",stds_w_area_values$eUnit),]
stds_w_area_values <- stds_w_area_values[-grep("RI.*",stds_w_area_values$eUnit),]

nrow(stds_w_area_values)
# [1] 1705

# checking how many different samples are covered by each of the standards
length(unique(stds_w_area_values[which(stds_w_area_values$Name == "2-hexanone-1,1,1,3.3-d5"), 'eUnit']$eUnit))
# [1] 586
length(unique(stds_w_area_values[which(stds_w_area_values$Name == "Benzaldehyde-d6"), 'eUnit']$eUnit))
# [1] 566
length(unique(stds_w_area_values[which(stds_w_area_values$Name == "Ethyl acetate-d8"), 'eUnit']$eUnit))
# [1] 553

# what is the distribution of area values
summary(as.numeric(stds_w_area_values$Area))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6103  947137 1136632 1746089 2847293 6079598

# checking to see how many duplicate pairs of nursery id and standard pairs are
# there.
cmpd_eunit <- paste(stds_w_area_values$Name, stds_w_area_values$eUnit)

stds_w_area_values <- cbind(stds_w_area_values, cmpd_eunit)

tail(sort(table(stds_w_area_values$cmpd_eunit)))
# there seem to be none, which is good.

# INTERPRETATION: This above table shows the the number of times any particular
# row is duplicated (meaning that a single sample spans multiple rows). We can
# see the rows are always singular, and there are no duplicates.

###################################################
## 2. GENERATING FINAL STANDARDIZATION DATAFRAME ##
###################################################

# getting the nursery id from eunit string
stds_w_area_values$nursery_id <- as.numeric(gsub(
  "eunit ([0-9]{1,})",
  "\\1",
  str_extract(stds_w_area_values$eUnit, "eunit [0-9]{1,}")
))

final_std.df <- stds_w_area_values[,c('Name','Area','nursery_id')]

# convert Area column to numeric
final_std.df$Area <- as.numeric(final_std.df$Area)

# keeping only the standards for which the Area values are over 10,000
final_std.df <- final_std.df[final_std.df$Area > 23000,]

# checking how many different samples are covered by each of the standards AGAIN
length(unique(final_std.df[which(final_std.df$Name == "2-hexanone-1,1,1,3.3-d5"),'nursery_id']))
# [1] 560
length(unique(final_std.df[which(final_std.df$Name == "Benzaldehyde-d6"),'nursery_id']))
# [1] 561
length(unique(final_std.df[which(final_std.df$Name == "Ethyl acetate-d8"),'nursery_id']))
# [1] 553


# generating a table of nursery_id and area values for the three standards
# something like this:
#
#  +------------+----------+-----------+------------+
#  | nursery_id | area_hex | area_benz | area_ethyl |
#  +------------+----------+-----------+------------+
#  |   ...      | ...      |  ...      | ...        |
#  +------------+----------+-----------+------------+
#

unique_nids <- unique(final_std.df$nursery_id)


nursery_id_std_area_vals <- NULL
for(nid in unique_nids) {
  area_val_hex <- final_std.df[
    which(final_std.df$nursery_id == nid &
            final_std.df$Name == "2-hexanone-1,1,1,3.3-d5"), 'Area']
  
  area_val_benz <- final_std.df[
    which(final_std.df$nursery_id == nid &
            final_std.df$Name == "Benzaldehyde-d6"), 'Area']
  
  area_val_ethyl <- final_std.df[
    which(final_std.df$nursery_id == nid &
            final_std.df$Name == "Ethyl acetate-d8"), 'Area']
  
  # set compound area to NA if the compound is not found in that nursery id
  if (length(area_val_hex) == 0) area_val_hex <- NA
  if (length(area_val_benz) == 0) area_val_benz <- NA
  if (length(area_val_ethyl) == 0) area_val_ethyl <- NA
  
  nursery_id_std_area_vals <- rbind(
    nursery_id_std_area_vals,
    cbind(nid, area_val_hex, area_val_benz, area_val_ethyl)
  )
}

# Checking the number of samples covered by all three standards
nrow(nursery_id_std_area_vals)
# [1] 561

nrow(
  as.data.frame(!is.na(nursery_id_std_area_vals)) %>% 
  filter(area_val_hex == TRUE, area_val_benz == TRUE, area_val_ethyl == TRUE)
)
# [1] 552

# INTERPRETATION: There are 552 samples that are covered by all the three
# standards.


# checking how many different missing area values are there for each of the 
# standards
sum(is.na(nursery_id_std_area_vals[,2]))
# [1] 1
sum(is.na(nursery_id_std_area_vals[,3]))
# [1] 0
sum(is.na(nursery_id_std_area_vals[,4]))
# [1] 8

summary(nursery_id_std_area_vals[,2])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 1609652 2872826 3295122 3312535 3697735 6079598       2

summary(nursery_id_std_area_vals[,3])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 523296  972429 1061802 1076549 1156829 1891278

summary(nursery_id_std_area_vals[,4])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  233084  765661  939015  936108 1108099 1774002       9

boxplot(
  nursery_id_std_area_vals[,2],
  nursery_id_std_area_vals[,3],
  nursery_id_std_area_vals[,4]
)

# Generate the distribution of Benzaldehyde standard
shapiro.test(nursery_id_std_area_vals[,3])
# Shapiro-Wilk normality test
# 
# data:  nursery_id_std_area_vals[, 3]
# W = 0.97427, p-value = 2.304e-08

plots <- list()

plots[[1]] <- ggqqplot(nursery_id_std_area_vals[,3]) + GLOBAL_THEME
plots[[2]] <- gghistogram(nursery_id_std_area_vals[,3]) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Volatile abundance (TIC)") + ylab("Count") +
  GLOBAL_THEME

annotate_figure(ggarrange(plotlist = plots,nrow=1,ncol=2))

nrow(nursery_id_std_area_vals)
# 561

# writing the nursery id standard area values table
colnames(nursery_id_std_area_vals) <- c(
  "nursery_id",
  "2-hexanone-1,1,1,3.3-d5",
  "Benzaldehyde-d6",
  "Ethyl acetate-d8"
)
 
write.table(
  nursery_id_std_area_vals,
  'data/processed/standardized/nursery_id_std_area_vals.tsv',
  sep = "\t"
)

########################################
## 3. GENERIC STANDARDIZATION ROUTINE ##
########################################

KNOWN_STANDARDS <- c("2-hexanone-1,1,1,3.3-d5", "Benzaldehyde-d6", "Ethyl acetate-d8")

# remove one of the replicates chosen at random
rand_rep_rem <- c(6167,1009,1197,2250,9304,9132,8087)
idx_rem <- NULL
for(i in seq_along(rand_rep_rem)){
  # gather the indexes to remove from standard df
  idx_rem <- c(idx_rem,which(final_std.df$nursery_id == rand_rep_rem[i]))
}

# rows before removing the replicates
nrow(final_std.df)
# [1] 1674

final_std.df <- final_std.df[-idx_rem,]

# rows after removing replicates
nrow(final_std.df)
# [1] 1647

standardize_with <- function(pheno_tbl, standard_name) {

  # check to see if the standard name is appropriate.
  if (standard_name %in% KNOWN_STANDARDS) {
    
    # moving nursery ids from colum to rownames
    nurseryids <- pheno_tbl$nurseryid
    pheno_tbl <- pheno_tbl[, !names(pheno_tbl) %in% "nurseryid"]
    rownames(pheno_tbl) <- nurseryids
    
    # get only the rows which belong to standard_name
    fnl_std_tbl <-
      final_std.df[
        which(
          final_std.df$Name == standard_name &
            !is.na(final_std.df$nursery_id)
        ),]
  
    nursery_ids <- unique(fnl_std_tbl$nursery_id)
    
    new_pheno_tbl <- matrix(NA, length(nursery_ids), ncol(pheno_tbl))
    rownames(new_pheno_tbl) <- nursery_ids
    colnames(new_pheno_tbl) <- colnames(pheno_tbl)
    
    area_value <- NULL
    # go through each nursery id and for that corresponding row in phenotypic
    # table and standardize all the columns of that row by the Area value.
    for (i in seq_along(nursery_ids)) {

      nid <- as.numeric(nursery_ids[i])

      area_value[i] <- as.numeric(fnl_std_tbl[which(nursery_ids == nid),2])

      
      new_pheno_tbl[i,] <-
        as.numeric(
          pheno_tbl[which(row.names(pheno_tbl) == nid),]
        ) / as.numeric(area_value[i])
    }
    
    new_pheno_tbl <- na.omit(new_pheno_tbl)
    
    return(new_pheno_tbl)
  } else {
    stop(
      standard_name, " is not a internal standard. Your options are: \n\t- ",
      paste(KNOWN_STANDARDS, collapse = "\n\t- ")
    )
  }

}

# import phenotype table
pheno_tbl <- read.table(
  'data/processed/phenotype_tables/gcms_phenotype_table_w_nurseryid_sim_gt_850.tsv',
  sep = "\t",
  header = TRUE
)

pheno_tbl_106 <- read.table(
  'data/processed/phenotype_tables/gcms_phenotype_table_w_nurseryid_manually_curated_106_cmpds.tsv',
  sep = "\t",
  header = TRUE
)

# remove one of the replicates chosen at random
rand_rep_rem <- c(6167,1009,1197,2250,9304,9132,8087)
idx_rem <- NULL
for(i in seq_along(rand_rep_rem)){
  # gather the indexes to remove from standard df
  idx_rem <- c(idx_rem,which(rownames(pheno_tbl_106) == rand_rep_rem[i] ))
}
# remove the replicates if they exist
if(length(idx_rem) > 0){
  # rows before removing the replicates
  nrow(pheno_tbl_106)
  # [1] 662
  
  pheno_tbl_106 <- pheno_tbl_106[-idx_rem,]

}
nursery_id <- rownames(pheno_tbl_106)
write.table(
  cbind(nursery_id,pheno_tbl_106),
  'data/processed/gcms_phenotype_table_w_nurseryid_no_replicates_sim_gt_850.tsv',
  sep = "\t",
  row.names = FALSE
)

#################################################
## 4. STANDARDIZE WITH 2-hexanone-1,1,1,3.3-d5 ##
#################################################
pheno_tbl_std_w_hexanone <-
  standardize_with(pheno_tbl, "2-hexanone-1,1,1,3.3-d5")

nurseryid <- rownames(pheno_tbl_std_w_hexanone)
pheno_tbl_std_w_hexanone <- cbind(nurseryid,pheno_tbl_std_w_hexanone)
write.table(
  pheno_tbl_std_w_hexanone,
  'data/processed/standardized/gcms_phenotype_table_std_w_hexanone_sim_gt_850.tsv',
  sep = "\t",
  row.names = FALSE
)

#########################################
## 5. STANDARDIZE WITH Benzaldehyde-d6 ##
#########################################
pheno_tbl_std_w_benzaldehyde <-
  standardize_with(pheno_tbl, "Benzaldehyde-d6")

nursery_id <- rownames(pheno_tbl_std_w_benzaldehyde)
pheno_tbl_std_w_benzaldehyde <- cbind(nursery_id,pheno_tbl_std_w_benzaldehyde)
write.table(
  pheno_tbl_std_w_benzaldehyde,
  'data/processed/standardized/gcms_phenotype_table_std_w_benzaldehyde_sim_gt_850.tsv',
  sep = "\t",
  row.names = FALSE
)

# standardizing the 109 compounds phenotype table
pheno_tbl_109_std_w_benzaldehyde <- 
  standardize_with(pheno_tbl_109, "Benzaldehyde-d6")

nurseryid <- rownames(pheno_tbl_109_std_w_benzaldehyde)
pheno_tbl_109_std_w_benzaldehyde <- cbind(nurseryid, pheno_tbl_109_std_w_benzaldehyde)
write.table(
  pheno_tbl_109_std_w_benzaldehyde,
  'data/processed/standardized/gcms_phenotype_table_std_w_benzaldehyde_sim_gt_850_109_cmpds.tsv',
  sep = "\t",
  row.names = FALSE
)

# standardizing the 106 compounds phenotype table
pheno_tbl_106_std_w_benzaldehyde <- 
  standardize_with(pheno_tbl_106, "Benzaldehyde-d6")

nurseryid <- rownames(pheno_tbl_106_std_w_benzaldehyde)
pheno_tbl_106_std_w_benzaldehyde <- cbind(nurseryid, pheno_tbl_106_std_w_benzaldehyde)
write.table(
  pheno_tbl_106_std_w_benzaldehyde,
  'data/processed/standardized/gcms_phenotype_table_std_w_benzaldehyde_sim_gt_850_106_cmpds.tsv',
  sep = "\t",
  row.names = FALSE
)


##########################################
## 6. STANDARDIZE WITH Ethyl acetate-d8 ##
##########################################
pheno_tbl_std_w_ethylacetate <-
  standardize_with(pheno_tbl, "Ethyl acetate-d8")

nursery_id <- rownames(pheno_tbl_std_w_ethylacetate)
pheno_tbl_std_w_ethylacetate <- cbind(nursery_id,pheno_tbl_std_w_ethylacetate)
write.table(
  pheno_tbl_std_w_ethylacetate,
  'data/processed/standardized/gcms_phenotype_table_std_w_ethylacetate_sim_gt_850.tsv',
  sep = "\t",
  row.names = FALSE
)

