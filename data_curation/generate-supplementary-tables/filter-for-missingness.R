# Title     : Filtered Supplemnetary Tables
# Objective : This script takes the generated supplementary tables and performs
#             missingness analysis, filters the data for certain missingness threshold
#             and exports the filtered data


############################
## IMPORTS & DATA LOADING ##
############################

library(readxl)
library(xlsx)
library(dplyr)

# load the GCMS dataset
gcms_data <- read_excel(
  'data/processed/sup_tbl_1-final_gcms_phenotype_table.xlsx'
)
# checking to see what the data structure is
dim(gcms_data)
# [1] 553 107
# There are 553 samples and 106 compounds + the first column which is the appleid

pheno_tbl <- read_excel(
  'data/processed/sup_tbl_2-abc_phenotype_table.xlsx'
)
dim(pheno_tbl)
# [1] 553 44
# Therea re 553 samples and 39 quantitative + 5 categorical phenotypes

# load the apple ids to keep that we have genetic data for
keep_aids <- read.table(
  'data/raw/appleids-with-genetic-data.txt'
)
keep_aids <- keep_aids$V1
length(keep_aids)
# [1] 516
# There are 516 apple ids that we would like to keep. The reset 37 are the ones that we do not have
# genetic data available.

# load the compound classification data
compound_classification_data <- read_excel(
  'data/processed/final_compounds_table_rev2_addressed.xlsx'
)

dim(compound_classification_data)
# [1] 106 5

# The above dimentions indicate that there are 106 compounds (as rows)
# and 5 columns.

colnames(compound_classification_data)
# [1] "Name"           "Samples (%)"    "Classification" "Likely (Y/N)"   "Comments"

############################################
## GENERATE THE KEEP & THROW AWAY DATASET ##
############################################

# extract the throw away apple ids
all_aids <- unlist(gcms_data$appleid)

# generate the keep dataset
keep_aids_idx <- which(gcms_data$appleid %in% keep_aids)
keep.gcms.df <- gcms_data[keep_aids_idx,]
keep.gcms.df.noaid <- keep.gcms.df[, 2:ncol(keep.gcms.df)]

# generate the throw away dataset
bad_aids <- setdiff(all_aids, keep_aids)
bad_aids_idx <- which(gcms_data$appleid %in% bad_aids)
throw_away.gcms.df <- gcms_data[bad_aids_idx,]
throw_away.gcms.df.noaid <- throw_away.gcms.df[, 2:ncol(throw_away.gcms.df)]

##########################################
## INVESTIGATING THE THROW AWAY DATASET ##
##########################################

dim(throw_away.gcms.df.noaid)
# [1]  37 106
# There are 37 accessions that we are throwing away

# figuring out the names of these apples that we do not have the genetic data for
bad_apples_idx <- which(pheno_tbl$apple_id %in% bad_aids)
throw_away.pheno.df <- pheno_tbl[bad_apples_idx,]

table(throw_away.pheno.df$origin)
# Kentville      USDA
# 31         6

table(throw_away.pheno.df$use)
# cider dessert    wild
# 1      28       2

# looking at the missingness of this throw away dataset

# missingness by samples
t_missingness_by_samples <- rowSums(throw_away.gcms.df.noaid == 0)
# [1] 54 53 44 60 40 65 59 56 54 50 56 69 66 39 63 45 67 53 46 51 51 62 68 56 51
# [26] 60 66 59 49 57 63 52 51 57 60 58 40
# This shows the number of compounds that are absent for each of the samples that we throw away.
# When we translate these to percentages, we get:

t_missingness_by_samples.percent <- (t_missingness_by_samples / ncol(throw_away.gcms.df.noaid)) * 100
round(t_missingness_by_samples.percent, 2)
# [1] 50.94 50.00 41.51 56.60 37.74 61.32 55.66 52.83 50.94 47.17 52.83 65.09
# [13] 62.26 36.79 59.43 42.45 63.21 50.00 43.40 48.11 48.11 58.49 64.15 52.83
# [25] 48.11 56.60 62.26 55.66 46.23 53.77 59.43 49.06 48.11 53.77 56.60 54.72
# [37] 37.74

# generating a histogram of the above distribution
hist(
  t_missingness_by_samples.percent, breaks = 30,
  main = "Distribution of missingness by samples",
  xlab = "Compounds Missing (%)",
  ylab = "Samples Count"
)

# missingness by compounds
t_missingness_by_compounds <- as.data.frame(colSums(throw_away.gcms.df.noaid == 0))
t_missingness_by_compounds$Compound <- rownames(t_missingness_by_compounds)
rownames(t_missingness_by_compounds) <- NULL
colnames(t_missingness_by_compounds) <- c("SamplesMissing", "Compound")
t_missingness_by_compounds$PercentSamplesMissing <- (t_missingness_by_compounds$SamplesMissing / 37) * 100

# generating the visualization of missingness by compounds
hist(
  t_missingness_by_compounds$PercentSamplesMissing, breaks = 50,
  main = "Distribution of missingness by compounds",
  xlab = "Samples Missing (%)",
  ylab = "Compounds Count"
)

# Q: In the above histogram, we see that there are 12 compound which have almost no missingness.
#    I am interested to know which compounds they are? Are the known to be significant in apples?
top_abundant_compounds <- t_missingness_by_compounds[which(t_missingness_by_compounds$PercentSamplesMissing == 0), 'Compound']
top_abundant_compounds
# [1] "1-Butanol"                  "1-Butanol, 2-methyl-"
# [3] "1-Hexanol"                  "1-Hexanol, 2-ethyl-"
# [5] "1-Penten-3-one"             "2-Heptenal"
# [7] "2-Hexenal"                  "2-4-Hexadienal"
# [9] "5-Hepten-2-one, 6-methyl-"  "Acetic-acid, hexyl-ester"
# [11] "Butanoic-acid, butyl-ester" "Hexanal"

# Here are these compounds. And of course, we have previously known that 1-Butanol nad 1-Hexanol are
# one of the most important compounds in apples and are quite abundant. This checks out and is a good
# evidence that we are on the right track.


# Q: I am further interested in knowing whether these compounds all belong to the
#    same classification of compounds or similar? Is there some relation between
#    those compounds? Are they all part of the same pathway?

top_abundant_compounds <- data.frame(Name = top_abundant_compounds)
top_abundant_compounds <- left_join(top_abundant_compounds, compound_classification_data, by = "Name")

table(top_abundant_compounds$Classification)
#
# Alcohol              Aldehyde Ester, Straight chain
# 4                     4                     2
# Ketone
# 2

# Okay, so there seem to be equal amounts of Alcohol and Esters. This is fitting
# because Esters are the most abundant VOCs in apples according to following citation:

# Sugimoto, N., Engelgau, P., Jones, A. D., Song, J., & Beaudry, R. (2021).
#   Citramalate synthase yields a biosynthetic pathway for isoleucine
#     and straight- and branched-chain ester formation in ripening apple fruit.
# Proceedings of the National Academy of Sciences, 118(3), e2009988118.
# https://doi.org/10.1073/pnas.2009988118

####################################
## INVESTIGATING THE KEEP DATASET ##
####################################

dim(keep.gcms.df.noaid)
# [1] 516 106
# There are 516 accessions that we are left with after throwing away the data for which we do not have the genetic data
# available

# missingness by samples
k_missingness_by_samples <- rowSums(keep.gcms.df.noaid == 0)
k_missingness_by_samples.percent <- (k_missingness_by_samples / ncol(keep.gcms.df.noaid)) * 100

hist(
  k_missingness_by_samples.percent, breaks = 100,
  main = "Distribution of missingness by samples",
  xlab = "Compound Missing (%)",
  ylab = "Sample Count"
)

# investigating the troubling sample which 90% missingness
troubling_sample_idx <- which(k_missingness_by_samples.percent >= 80)
troubling_aid <- as.numeric(keep.gcms.df[troubling_sample_idx, 'appleid'])
pheno_tbl[which(pheno_tbl$apple_id == troubling_aid), c('world', 'use', 'origin', 'country')]
# world use   origin country
# <chr> <chr> <chr>  <chr>
# 1 new   dessert USDA   United States

# looking at the missingness of the ABC phenotypes
pheno_idxs <- 6:ncol(pheno_tbl)
num_phenos_missing <- sum(is.na(pheno_tbl[which(pheno_tbl$apple_id == troubling_aid), pheno_idxs]))
pcent_phenos_missing <- (num_phenos_missing / length(pheno_idxs)) * 100
pcent_phenos_missing
# [1] 66.6

# This outlier is a German crab apple most of whose ABC phenotypes are also missing. Therefore
# it is safe to get rid of this.

# removing the troubling outlier sample
keep.gcms.df <- keep.gcms.df[-troubling_sample_idx,]
keep.gcms.df.noaid <- keep.gcms.df[, 2:ncol(keep.gcms.df)]

dim(keep.gcms.df.noaid)
# [1] 515 106

# we need to generate the missingness by samples histogram again to see what the distribution looks like.
k_missingness_by_samples <- rowSums(keep.gcms.df.noaid == 0)
k_missingness_by_samples.percent <- (k_missingness_by_samples / ncol(keep.gcms.df.noaid)) * 100

hist(
  k_missingness_by_samples.percent, breaks = 100,
  main = "Distribution of missingness by samples",
  xlab = "Compound Missing (%)",
  ylab = "Sample Count"
)

# missingness by compounds
k_missingness_by_compounds <- as.data.frame(colSums(keep.gcms.df.noaid == 0))
k_missingness_by_compounds$Compound <- rownames(k_missingness_by_compounds)
rownames(k_missingness_by_compounds) <- NULL
colnames(k_missingness_by_compounds) <- c("SamplesMissing", "Compound")
k_missingness_by_compounds$PercentSamplesMissing <-
  (k_missingness_by_compounds$SamplesMissing / nrow(keep.gcms.df.noaid)) * 100

# generating the visualization of missingness by compounds
hist(
  k_missingness_by_compounds[which(k_missingness_by_compounds$PercentSamplesMissing >= 4), 'PercentSamplesMissing'], breaks = 50,
  main = "Distribution of missingness by compounds",
  xlab = "Samples Missing (%)",
  ylab = "Compounds Count"
)
summary(k_missingness_by_compounds$PercentSamplesMissing)
# summary(k_missingness_by_compounds$PercentSamplesMissing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00   20.49   70.00   55.51   84.81   93.01

k_missingness_by_compounds[which(k_missingness_by_compounds$PercentSamplesMissing < 4), 'Compound']
# [1] "1-Butanol"                  "1-Butanol, 2-methyl-"
# [3] "1-Hexanol"                  "1-Hexanol, 2-ethyl-"
# [5] "2-Hexenal"                  "2-4-Hexadienal"
# [7] "5-Hepten-2-one, 6-methyl-"  "a-Farnesene"
# [9] "Acetic-acid, butyl-ester"   "Acetic-acid, hexyl-ester"
# [11] "Butanoic-acid, butyl-ester" "Hexanal"

# In conclusion, the final dimensions for the GCMS phenotype dataset are:
dim(keep.gcms.df)
# [1] 515 107
# NOTE: Extra column is for apple id.

# export the absolutely final GC-MS phenotype dataset.
openxlsx::write.xlsx(
  keep.gcms.df, "data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx"
)

# export the version 2 of ABC phenotype table that only contains the same apples as the
# one in the GCMS data
abc_pheno_tbl <- read_excel(
  'data/processed/sup_tbl_2-abc_phenotype_table.xlsx'
)

final_keep_aids <- keep.gcms.df$appleid
abc_keep_aids_idx <- which(abc_pheno_tbl$apple_id %in% final_keep_aids)

# new ABC phenotype data
abc_pheno_tbl <- abc_pheno_tbl[abc_keep_aids_idx,]

openxlsx::write.xlsx(
  abc_pheno_tbl, 'data/processed/sup_tbl_2-abc_phenotype_table_v2.xlsx'
)
