
library(hash)
library(tidyverse)
library(data.table)
library(zeallot)
library(ape)

source('exploration/gwas_summary/utils.R')

##################
## DATA LOADING ##
##################

class_corrected_gwas_compound_files <- list.files(
    '/Users/tayabsoomro/Documents/Projects/GCMS-Project/class_corrected_gwas_analysis/gwas_data_updated_names/all-samples/',
    pattern = "*/*std_pvals.csv", recursive = TRUE, full.names = TRUE
)
length(class_corrected_gwas_compound_files)
# [1] 101

gwas_1_compound_files <- list.files(
  '/Users/tayabsoomro/Documents/Projects/GCMS-Project/gwas_1/gwas_results/',
  pattern = "*/*_std_pvals.csv", recursive = TRUE, full.names = TRUE
)
length(gwas_1_compound_files)
# [1] 104

# load the GFF3 file
gff3_dat <- read.gff('exploration/gwas_summary/Malus_x_domestica.v1.0.consensus.gff')
dim(gff3_dat)
# [1] 508830      9

# only keep the mRNA annotations
gff3_dat <- gff3_dat[which(gff3_dat$type == "mRNA"),]
dim(gff3_dat)
# [1] 63541     9


###############################
## GATHER THE COMPOUND NAMES ##
###############################

compounds_to_compare_between_two_gwas <- c(
  "1-Butanol", "1-Hexanol", "2-methylbutyl acetate",
  "Butyl acetate", "Hexyl acetate", "Pentyl acetate"
)
length(compounds_to_compare_between_two_gwas)
# [1] 6

compounds_only_in_class_corrected_gwas <- c(
  "1-Pentanol", "Propyl butyrate"
)
length(compounds_only_in_class_corrected_gwas)
# [1] 2

compounds_only_in_gwas_1 <- c(
  "2-4-Hexadienal", "Isobutyl acetate", "n-Propyl-acetate"
)
length(compounds_only_in_gwas_1)
# [1] 3

all_compounds <- c(
  compounds_to_compare_between_two_gwas,
  compounds_only_in_class_corrected_gwas,
  compounds_only_in_gwas_1
)
length(all_compounds)
# [1] 11

#############################
## GATHER ALL THE TOP HITS ##
#############################

compound_map_peak <- hash()
top_hits <- data.frame()

# get the top hits for class_corrected gwas
c(top_hits,compound_map_peak) %<-% update_top_hits(
  compounds_to_compare_between_two_gwas,
  class_corrected_gwas_compound_files,
  top_hits, compound_map_peak
)
dim(top_hits)
# [1] 6 6

# get the top hits for class_corrected gwas_1
c(top_hits,compound_map_peak) %<-% update_top_hits(
  compounds_to_compare_between_two_gwas, 
  gwas_1_compound_files, 
  top_hits, compound_map_peak
)
dim(top_hits)
# [1] 12 6

# get the top hits for only class_corrected gwas
c(top_hits,compound_map_peak) %<-% update_top_hits(
  compounds_only_in_gwas_1, 
  gwas_1_compound_files, 
  top_hits, compound_map_peak
)
dim(top_hits)
# [1] 15 6

# get the top hits for only gwas_1
c(top_hits,compound_map_peak) %<-% update_top_hits(
  compounds_only_in_class_corrected_gwas, 
  class_corrected_gwas_compound_files, 
  top_hits, compound_map_peak
)
dim(top_hits)
# [1] 17 6


###############################################################
## GATHER ANNOTATION FOR ANYTHING +/- 10kb FROM THE TOP HITS ##
###############################################################

# add the attributes to the top_hits data frame
attributes_list <- c()
for(i in seq_len(nrow(top_hits))){
  # get the lower and upper bound of the position
  pos <- as.numeric(top_hits[i,"POS"])
  chr <- paste0("chr",as.numeric(top_hits[i,"CHR"]))
  print(chr)
  window_left <- pos - 100000
  window_right <- pos + 100000
  attributes <- gff3_dat[
    which(gff3_dat$seqid == chr & gff3_dat$start >= window_left & gff3_dat$end <= window_right),
    'attributes'
  ]
  attributes_str <- ""
  for(a in attributes){
    attributes_str <- paste0(attributes_str,";",a)
  }
  attributes_list <- c(attributes_list, attributes_str)
}

top_hits$attributes <- attributes_list

write.xlsx(
  top_hits,
  "exploration/gwas_summary/top_hits.xlsx"
)