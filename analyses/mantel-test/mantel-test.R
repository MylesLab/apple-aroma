# Title     : Mantel Test
# Objective : This script performs correlations between genetic and phenotype
#             data
# Created by: tayabsoomro
# Created on: 2021-06-05

############################
## IMPORTS & DATA LOADING ##
############################

library(dplyr)
library(ade4)
library(readxl)

gcms_pheno_tbl <- read_excel(
  'data/processed/Supplementary_Data.xlsx',
  sheet = 'GCMS Data'
)
dim(gcms_pheno_tbl)
# [1] 515 107

# loading the genetic distance matrix
genetic_dist_mat <- read.table(
  'data/processed/ibs/plink.mdist'
)

# import the apple IDs from PLINK
ids <- read.table(
  'data/processed/ibs/ids.txt'
)
colnames(ids) <- "appleid"

# update the colnames and rownames for the genetic data matrix
colnames(genetic_dist_mat) <- as.numeric(ids$appleid)
rownames(genetic_dist_mat) <- as.numeric(ids$appleid)

##########################################
## PREPARE THE GENOTYPE DISTANCE MATRIX ##
##########################################

# keep only the samples that are present in our phenotype table
keep_aids <- unlist(gcms_pheno_tbl$appleid)

genetic_dist_mat <- genetic_dist_mat[
  rownames(genetic_dist_mat) %in% keep_aids,
  colnames(genetic_dist_mat) %in% keep_aids
]

###########################################
## PREPARE THE PHENOTYPE DISTANCE MATRIX ##
###########################################

# generate the distance matrix for phenotype table
rownames(gcms_pheno_tbl) <- unlist(gcms_pheno_tbl$appleid)
gcms_pheno_tbl.noaid <- gcms_pheno_tbl[, 2:ncol(gcms_pheno_tbl)]
gcms_dist_mat <- as.data.frame(as.matrix(dist(gcms_pheno_tbl.noaid)))

rownames(gcms_dist_mat) <- as.character(gcms_pheno_tbl$appleid)
colnames(gcms_dist_mat) <- as.character(gcms_pheno_tbl$appleid)

# making sure that the columns and rows in phenotype distance matrix align with
# the columns and rows in genetic distance matrix

gcms_dist_mat <- gcms_dist_mat[
  match(rownames(genetic_dist_mat), rownames(gcms_dist_mat)),
  match(colnames(genetic_dist_mat), colnames(gcms_dist_mat))
]

mantel.rtest(as.dist(genetic_dist_mat), as.dist(gcms_dist_mat), nrepet = 10000)
# Monte-Carlo test
# Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
#
# Observation: 0.05998306
#
# Based on 10000 replicates
# Simulated p-value: 0.00789921
# Alternative hypothesis: greater
#
# Std.Obs  Expectation     Variance
# 2.5531952272 0.0001401793 0.0005493599

