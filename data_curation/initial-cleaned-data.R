# Title     : Initial Data Cleanup
# Objective : Cleans up the initial data received

#######################
## TABLE OF CONTENTS ##
#######################

## 0. IMPORTS & DATA LOADING
## 1. REMOVING THE SATURATED ENCLOSING
## 2. REMOVING BLANKS
## 3. REMOVING DUTERATED COMPOUNDS
## 4. REMOVING ARTIFACTS
## 5. APPENDING APPLE ID COLUMN
## 6. REMOVING COMPOUNDS BELOW CERTAIN SIMILARITY THRESHOLD


library(dplyr)
library(readxl)
library(stringr)

###############################
## 0. IMPORTS & DATA LOADING ##
###############################

# Filename Variables
gcms_file <- 'data/raw/gcms_data.xlsx'
aid_nid_pivot_file <- 'data/raw/nursery-id_apple-id_pivot.tsv'

# Data Frames
gcms.df <- read_excel(gcms_file, col_types = "text")
dim(gcms.df)
# [1] 130376     12

pivot.df <- read.table(aid_nid_pivot_file)
dim(pivot.df)
# [1] 2356    2

#########################################
## 1. REMOVING THE SATURATED ENCLOSING ##
#########################################

# We realize that in the Area column, some of the values are enclosed by
# a saturated keyword. This happens because the GC-MS machine detects samples
# that are so highly abundant, that they saturate the entire column. Therefore,
# the values enclosing the saturated represent the highest possible Area value
# before the machine could not quantify anymore. We would like to keep those
# values because they convey how abundant that particular compound is. As such, we
# would like to peel out the "saturated(" and ")" part from the Area column
# if one exists.

# First, we are checking to see what the class type is in the Area value. We
# expect this to be "character" because of the "saturated" keyword.
class(gcms.df$Area)
# [1] "character"
# This ensures that the columns (particularly Area column) in the gcms.df are
# represented as character. This is important for removing the "saturated"
# keyword.

sum(is.na(gcms.df$Area))
# [1] 1865
# There are 1,865 rows that are not available in Area column.

unique(gcms.df[which(is.na(gcms.df$Area)), 'Area'])
# A tibble: 1 × 1
# Area
# <chr>
# 1 NA
# All of those above rows have the value "NA" in the Area column

# getting all the rows for which the Area field is saturated
gcms.df[grep("sat.*", gcms.df$Area),]$Area
# There are 160 rows which have saturated Area values

# Removing the saturated enclosing from the numbers in the Area field
sat_removed <- gsub('saturated\\( ', '', gcms.df$Area)
finally_sat_removed <- gsub('\\)', '', sat_removed)
gcms.df$Area <- as.numeric(finally_sat_removed)


# check that there are no saturated enclosings anymore.
gcms.df[grep("sat.*", gcms.df$Area),]$Area
# numeric(0)
# Confirmed that all saturated enclosings are removed.


########################
## 2. REMOVING BLANKS ##
########################

# It was verified by our AAFC collaborators that the only rows in the gcms_df
# that we need to keep are the ones that have the pattern 'eunit
# {Number}-{Number}' in the file name column. The number folloiwng the
# keyword "eunit" belongs to the nursery IDs for which the particular sample
# belonged to. Everything else that remains are the blanks and we are not
# interested in those for our analysis.

dim(gcms.df)
# [1] 130376     12

# getting all the rows that have eunit keyword in the filename column
gcms.df <- gcms.df[grep("eunit", gcms.df$FileName),]
dim(gcms.df)
# [1] 124850     12

# After keeping only the samples that are not blanks (i.e., have eunit keyword)
# there are still other samples that need to be removed. These are the rows
# for which the compound name is "Unknown". We are not interested in these as
# they cannot tell us any useful information.

# checking to see how many Unknowns are there
length(grep("Unknown", gcms.df$Name))
# [1] 12484

# removing all the rows for which the compound name is Unknown
gcms.df <- gcms.df[grep("Unknown", gcms.df$Name, invert = TRUE),]
dim(gcms.df)
# [1] 112366     12

#####################################
## 3. REMOVING DUTERATED COMPOUNDS ##
#####################################

# Duteration is a process whereby a hydrogen is replaced with deuterium. Duterated
# compounds are used in GCMS as internal standards. These compounds need to be
# removed from our analysis as they do not contribute to the analysis. These
# compound names have pattern: *d[0-9].* (they contain anywhere in the compound
# name the letter d, followed any amount of numbers).

# getting the row index of all the duteration compounds in the gcms dataset.
duteration <- grep("*d[0-9].*", gcms.df$Name)

# obtaining a tally of all the duteration compounds and the libraries that they
# belong to.
table(gcms.df[duteration, c('Name', 'Library')])
# Name                      AAFC_library_2018 mainlib
# 2-hexanone-1,1,1,3.3-d5               763       0
# 2-phenyl-d5-ethanol                   563       0
# Acenaphthene-d10                        0     155
# Benzaldehyde-d6                       575       0
# Ethyl acetate-d8                      553       0

# The above compounds were checked with our AAFC collaborators to ensure that
# they are in fact the ones to be removed.

# removing all the duteration compounds from gcms dataset
gcms.df <- gcms.df[-duteration,]

# ensuring that the compounds were successfully removed.
duterated_compound_names <- c(
  "2-hexanone-1,1,1,3.3-d5",
  "2-phenyl-d5-ethanol",
  "Acenaphthene-d10",
  "Benzaldehyde-d6", "Ethyl acetate-d8"
)

# looping through all duterated compounds to make sure that each of them
# should not be found in the Name column of gcms.df
found <- NULL
for (cmpd in duterated_compound_names) {
  found <- c(found, gcms.df$Name == cmpd)
}
unique(found)
# [1] FALSE    NA

# All of them are FALSE, hence the duteration removal process was succesful.

############################
## 4. REMOVING ARTIFACTS ##
############################

# Reading in list of artifact compounds and removing them from the
# gcms.df

# loading the data containing artifacts to remove
artifacts <- read_excel('./data/raw/artifacts_to_remove.xlsx')
artifacts.names <- artifacts$Name

length(artifacts.names)
# [1] 26

# cleaning the extra characters that are around the names
artifacts.names <- gsub(
  ".*\\] ",
  "",
  gsub(
    "\"(.*)\"",
    "\\1",
    artifacts.names
  )
)

# checking how many occurrences of the artifact names are there in the GCMS dataset.
num_occur <- NULL
for (i in seq_along(artifacts.names)) {
  num_occur <- c(num_occur, sum(gcms.df$Name %in% artifacts.names[i]))
}

num_occur
# [1] 12103   753   246  1089   603    55   402   112    57     3     5    23
# [13]    19    37    20     1     4     6     3     1     1     1     1     1
# [25]     1     1
# We can see that all of the artifacts in the file are at least present once.

# this above is just a sanity check to ensure that each artifact is
# actually present at least once. If there were zeros, that would mean that some
# of the compounds were not matched with the data.

length(which(gcms.df$Name %in% artifacts.names))
# [1] 15548
# there are 15,548 rows in which there are one of those 26 compounds that are
# known as artifacts and should be removed. NOTE that the reason why we have
# way more entries than 26 is because those 26 artifacts are unique within the
# initial dataset

# checking the number of rows before removing the artifacts
nrow(gcms.df)
# [1] 109,757

# checking the number of rows that will be removed
nrow(gcms.df[!(gcms.df$Name %in% artifacts.names),])
# [1] 94,209

# removing all the artifacts from the gcms dataset.
gcms.df <- gcms.df[!(gcms.df$Name %in% artifacts.names),]
dim(gcms.df)
# [1] 94209    12

# The number of artifacts were 15548. 109757-94209 = 15548. Math works!

write.table(gcms.df, 'data/processed/cleaned/cleaned-gcms-data.tsv', sep = "\t")

##################################
## 5. APPENDING APPLE ID COLUMN ##
##################################

# We need to add apple id column into our GC-MS dataset so that further analyses
# can be carried out. This is going to be done through a pivot table which
# maps the apple ids with the nursery ids. The nursery id in the GCMS dataset
# takes the format of "eunit {Number}_{Number} in the FileName column.

# parsing out the nursery ids
nurseryIds <- gsub(
  'eunit ',
  '',
  str_extract(gcms.df$FileName, 'eunit [0-9]{1,}')
)

length(nurseryIds)
# [1] 94209

nrow(gcms.df)
# [1] 94209

# The nurseryIds vector is just as long as the number of rows in the gcms.df
# therefore, the both are ordered and can be attached together.

gcms.df$nursery_id <- as.numeric(nurseryIds)

# adding the apple_id column through the pivot table
final.df <- left_join(gcms.df, pivot.df, by = 'nursery_id')

# checking to ensure that the join worked correctly
nrow(gcms.df)
# [1] 94209

nrow(final.df)
# [1] 94209

# the number of rows is maintained which is what we expect to happen because
# it is joining on top of the GCMS dataset. Manually checked a few entries to
# see and make sure that the merging was correct.

# exporting the data into a data file
write.table(
  final.df, 'data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping.tsv',
  sep = "\t"
)

##############################################################
## 6. REMOVING COMPOUNDS BELOW CERTAIN SIMILARITY THRESHOLD ##
##############################################################

# In order to obtain the compound names after the GCMS analysis. The peaks are
# matched to National Institute of Standards and Technology (NIST) Mass
# Spectrometry Database. The match score that is obtained from this matching
# is referred to as the similarity score. It was advised by our AAFC
# collaborators that we discard every sample that has similarity score of
# less than 700.

sum(final.df$Similarity > 800)
# [1] 57519

sum(final.df$Similarity > 600)
# [1] 94083

# checking the histogram for the similarity
hist(
  as.numeric(final.df$Similarity),
  xlab = "Similarity Score", main = "Compound Density based on similarity score"
)
# We noticed that there are not a lot of compounds for which the similarity is
# less than 600, however there are quite a few compounds where the similarity
# is less than 700. From the histogram, it looks like the compounds have been
# somewhat cleaned to only include similarity of almost greater than 600.

# exporting the data into a data file
save_cleaned_dat_sim_thres <- function(sim_score) {
  write.table(
    final.df[which(as.numeric(final.df$Similarity) > sim_score),],
    paste0('data/processed/cleaned/sim_threshold/cleaned-gcms-data-with-appleid-mapping-sim-gt-',sim_score,'.tsv'),
    sep = "\t",
    row.names = FALSE
  )
}

sim_thres_buckets <- seq(600,950,50)
for (score in sim_thres_buckets) {
  save_cleaned_dat_sim_thres(score)
}


