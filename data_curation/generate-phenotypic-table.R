# Title     : Generate Phenotypic Table
# Objective : This script is responsible for generating the phenotypic table
#             GC-MS data.

# LIBRARIES
library(foreach)
library(dplyr)
library(doParallel)
library(MASS)

generate_pheno_tbl <- function(cleaned_data, w_aid_fname, w_nid_fname) {

  compoundNames <- unique(cleaned_data$Name) # 2875 compounds
  length(compoundNames)
  # [1] 2875
  
  appleIDs <- unique(cleaned_data$apple_id) # 663 apples
  length(appleIDs)
  # [1] 663
  
  nurseryIDs <- unique(cleaned_data$nursery_id) # 670 apples
  length(nurseryIDs)
  # [1] 670
  
  # NOTE: the reason why there is a difference of 7 between length of unique
  # apple ids and unique nursery ids is that there are 7 apple ids that belong
  # to multiple nursery ids.
  
  # cleaning up and only keeping relevant columns
  final_cleaned.df <- cleaned_data[, c('Name', 'apple_id', 'nursery_id', 'Area')]
  
  # check the columns
  colnames(final_cleaned.df)
  # [1] "Name"       "apple_id"   "nursery_id" "Area"
  
  # Combining all the replicates by taking the mean of their compound Area values.
  final_cleaned.df <- final_cleaned.df %>%
    group_by(Name, apple_id) %>%
    summarize(Area = mean(Area), nursery_id = min(nursery_id))
  
  # Now there are some inconsistencies between the nursery ids because some of
  # the rows will have nursery id from first replicate, and others will have it
  # from the second replicate. To clean that up, I am going to be using the
  # minimum values for the replicate nursery ids.
  
  bio_rep_nids <- list(
    c(6166, 6167),
    c(1008, 1009),
    c(1196, 1197),
    c(2249, 2250),
    c(9304, 9305),
    c(9132, 9134),
    c(8087, 8089)
  )
  
  # make sure that only the first nursery ids are used for all of the rows
  new_nursery_ids <- final_cleaned.df$nursery_id
  for (nid in bio_rep_nids) {
      new_nursery_ids <- gsub(nid[2], nid[1], new_nursery_ids)
  }
  final_cleaned.df$nursery_id <- as.numeric(new_nursery_ids)
  
  
  compoundNames <- unique(final_cleaned.df$Name) # 2875 compounds
  length(compoundNames)
  # [1] 2875
  
  appleIDs <- unique(final_cleaned.df$apple_id) # 663 apples
  length(appleIDs)
  # [1] 663
  
  nurseryIDs <- unique(final_cleaned.df$nursery_id) # 670 apples
  length(nurseryIDs)
  # [1] 662
  
  # the final tables to accumulate the data to
  pheno_tbl_aid <- matrix(, nrow = length(appleIDs), ncol = length(compoundNames))
  pheno_tbl_nid <- matrix(, nrow = length(nurseryIDs), ncol = length(compoundNames))
  
  # setting up parallel processing
  ncores <- detectCores()
  clust <- makeCluster(ncores - 1)
  registerDoParallel(clust)
  
  # accumulating the phenotype table with apple ids as rows
  pheno_tbl_aid <- foreach(aid = appleIDs, .combine = 'rbind') %:%
    foreach(cname = compoundNames, .combine = 'c') %dopar% {
      index <- which(final_cleaned.df$apple_id == aid & final_cleaned.df$Name == cname)
  
      sum(final_cleaned.df[index, 'Area'])
  }
  
  # add the column and row names
  colnames(pheno_tbl_aid) <- compoundNames
  rownames(pheno_tbl_aid) <- appleIDs
  
  # accumulating the phenotype table with nursery ids as rows
  pheno_tbl_nid <- foreach(nid = nurseryIDs, .combine = 'rbind') %:%
    foreach(cname = compoundNames, .combine = 'c') %dopar% {
      index <- which(final_cleaned.df$nursery_id == nid & final_cleaned.df$Name == cname)
  
      sum(final_cleaned.df[index, 'Area'])
  }
  
  # stop parallel processing
  stopCluster(clust)
  
  # add the column and row names
  colnames(pheno_tbl_nid) <- compoundNames
  row.names(pheno_tbl_nid) <- nurseryIDs
  
  pheno_tbl_nid <- cbind(nurseryid = nurseryIDs, pheno_tbl_nid)
  pheno_tbl_aid <- cbind(appleid = appleIDs, pheno_tbl_aid)
  
  # Writing the above tables to file
  write.table(
    pheno_tbl_aid, file = paste0('data/processed/phenotype_tables/',w_aid_fname, '.tsv'),
    sep = "\t",
    quote = TRUE,
    row.names = FALSE
  )
  
  write.table(
    pheno_tbl_nid, file = paste0('data/processed/phenotype_tables/',w_nid_fname, '.tsv'),
    sep = "\t",
    quote = TRUE,
    row.names = FALSE
  ) 
}

# Seeing how the phenotype table is affected with different similarity
# cutoffs.

# generate the phenotype tables for a similarity score threshold
generate_all_pheno_tbl_sim_score <- function(sim_score) {
  
  final_cleaned.df <- read.table(
    paste0(
      'data/processed/cleaned/sim_threshold/cleaned-gcms-data-with-appleid-mapping-sim-gt-',
      sim_score,'.tsv'
    ),
    header = TRUE
  )
  
  generate_pheno_tbl(
    final_cleaned.df,
    paste0('gcms_phenotype_table_w_appleid_sim_gt_',sim_score,'.tsv'),
    paste0('gcms_phenotype_table_w_nurseryid_sim_gt_',sim_score,'.tsv')
  )
}

# run the command to generate the phenotype files
sim_score_buckets <- seq(600,950,50)
for (score in sim_score_buckets) {
  print(score)
  generate_all_pheno_tbl_sim_score(score)
}

# generating the semi-final phenotype table after manually curated compounds
curated_cleaned.df <- read.table(
  'data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated_109.tsv',
  header = TRUE
)
generate_pheno_tbl(
  curated_cleaned.df,
  'gcms_phenotype_table_w_appleid_manually_curated_109_cmpds.tsv',
  'gcms_phenotype_table_w_nurseryid_manually_curated_109_cmpds.tsv'
)

# generating the final phenotype table after manually curated compounds and verified
# by AAFC collaborators.
final_curated_cleaned.df <- read.table(
  'data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated_106.tsv',
  header = TRUE
)

generate_pheno_tbl(
  final_curated_cleaned.df,
  'gcms_phenotype_table_w_appleid_manually_curated_106_cmpds',
  'gcms_phenotype_table_w_nurseryid_manually_curated_106_cmpds'
)
