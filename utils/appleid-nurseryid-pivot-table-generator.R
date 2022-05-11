# Title     : Apple ID & Nursery ID Pivot Table Generator.
# Objective : This script is responsible for generating the nursery id and
# pivot id table

# load the meta table of ABC
appleinfo.df <- read.table('data/raw/ABC_AppleInfo.csv', sep = ",", header = TRUE)
dim(appleinfo.df)
# [1] 2356   18

# see what columns are available
colnames(appleinfo.df)
# [1] "Row"             "Block"           "Tree"            "nursery_id"
# [5] "year_planted"    "guard_tree"      "check_tree"      "dead_tree"
# [9] "PLANTID"         "PLANTID_ALT1"    "ACCID"           "ACP"
# [13] "ACNO"            "ACS"             "kentville_block" "origin"
# [17] "apple_id"        "type"

# We need apple_id and nursery_id columns side by side
aid_nid_tbl <- appleinfo.df[, c('apple_id', 'nursery_id')]
# [1] 2356    2

# check to see if there are any duplicated rows
sum(duplicated(aid_nid_tbl))
# [1] 0

# there are no duplicated rows.

dim(aid_nid_tbl)
# [1] 2356    2

# exporting this table as our pivot table for nursery id and apple id
write.table(aid_nid_tbl, 'data/raw/nursery-id_apple-id_pivot.tsv', sep = "\t")
