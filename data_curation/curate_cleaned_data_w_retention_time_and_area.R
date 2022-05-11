# Title     : Retention Time & Area Values Based Cleanup
# Objective : This script generates a final dataset where the compound names
#             are curated and the nonsensical rows are removed based on the 
#             retention plots in figures/retention
# Created by: tayabsoomro
# Created on: 2021-05-20

library(dplyr)

# BACKGROUND: We generated retention plots for time vs time and area vs time to
# see whether the retention time is consistent for all the compounds. And we found
# that bulk of compounds had very similar retention time but there were some that
# were outliers. We need to remove those outliers. There were also instances 
# where there were two clusters of retention times. In such cases, we will consider
# them different compounds and the names of these should be changed with the 
# suffix of the first retention time value.

# Based on the retention plots where there were 16 plots in one image. We see 
# that after about 8 images we see that the number of samples decreases significantly
# and thus will not be usable for our downstream anlayses. Therefore, we will
# only curate the compounds on the first 8 images, which will yield 16*8=128 
# compounds

# import the cleaned dataset
cleaned.df <- read.table(
  'data/processed/cleaned/sim_threshold/cleaned-gcms-data-with-appleid-mapping-sim-gt-850.tsv',
  header = TRUE
)

###################
# 2,4-Hexadienal ##
###################

cleaned.df$Name <- gsub("2,4-Hexadienal, \\(E,E\\)-","2,4-Hexadienal",cleaned.df$Name)
hist(cleaned.df[which(cleaned.df$Name == "2,4-Hexadienal"),'X1st.Dimension.Time..s.'], breaks = 100)
idx_to_rem <- which(
  cleaned.df$Name == "2,4-Hexadienal" & 
  (cleaned.df$X1st.Dimension.Time..s. < 250 | cleaned.df$X1st.Dimension.Time..s. > 500) | 
    (cleaned.df$X1st.Dimension.Time..s. > 400 & cleaned.df$X1st.Dimension.Time..s. < 425)
)
cleaned.df <- cleaned.df[-c(idx_to_rem),]

# split 1 for this compound
hist(cleaned.df[which(cleaned.df$Name == "2,4-Hexadienal"),'X1st.Dimension.Time..s.'], breaks = 100)
first_mean <- mean(cleaned.df[which(cleaned.df$Name == '2,4-Hexadienal' & cleaned.df$X1st.Dimension.Time..s. < 295),'X1st.Dimension.Time..s.'])
cleaned.df[which(
  cleaned.df$Name == '2,4-Hexadienal' & cleaned.df$X1st.Dimension.Time..s. < 295
),'Name'] = '2,4-Hexadienal-294'

# split 2 for this compound
hist(cleaned.df[which(cleaned.df$Name == '2,4-Hexadienal'),'X1st.Dimension.Time..s.'])
second_mean <- mean(cleaned.df[which(cleaned.df$Name == '2,4-Hexadienal'),'X1st.Dimension.Time..s.'])
cleaned.df[which(
  cleaned.df$Name == '2,4-Hexadienal' & cleaned.df$X1st.Dimension.Time..s. > 400
),'Name'] = '2,4-Hexadienal-430'

# remove any residual compounds
cleaned.df <- cleaned.df[-which(cleaned.df$Name == "2,4-Hexadienal"),]

##############
# 2-Hexenal ##
##############

# change the name to be consistent through the compounds
cleaned.df[grep("2-Hexenal.*",cleaned.df$Name),'Name'] = '2-Hexenal'

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '2-Hexenal' & 
    (cleaned.df$X1st.Dimension.Time..s. < 293 | cleaned.df$X1st.Dimension.Time..s. > 295)
),]

##############
# 3-Hexenal ##
##############
cleaned.df[grep("3-Hexenal.*",cleaned.df$Name),'Name'] <- '3-Hexenal'

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '3-Hexenal' & 
    (cleaned.df$X1st.Dimension.Time..s. < 241 | cleaned.df$X1st.Dimension.Time..s. > 242)
),]

#############################
# Acetic acid, butyl ester ##
#############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Acetic acid, butyl ester' & 
    (cleaned.df$X1st.Dimension.Time..s. < 195 | cleaned.df$X1st.Dimension.Time..s. > 200)
),]

#############################
# Acetic acid, hexyl ester ##
#############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Acetic acid, hexyl ester' & 
    (cleaned.df$X1st.Dimension.Time..s. < 300 | cleaned.df$X1st.Dimension.Time..s. > 335)
),]

############
# Hexanal ##
############

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Hexanal' & 
    (cleaned.df$X1st.Dimension.Time..s. < 202 | cleaned.df$X1st.Dimension.Time..s. > 203)
),]

##############
# 1-Butanol ##
##############

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '1-Butanol' & 
    (cleaned.df$X1st.Dimension.Time..s. < 235 | cleaned.df$X1st.Dimension.Time..s. > 236)
),]

#########################
# 1-Butanol, 2-methyl- ##
#########################

# fix the name
cleaned.df[grep("1-Butanol, 2-methyl-, \\(S\\)-",cleaned.df$Name),'Name'] <- '1-Butanol, 2-methyl-'

###############################################################
# 2-Buten-1-one, 1-(2,6,6-trimethyl-1,3-cyclohexadien-1-yl)- ##
###############################################################

cleaned.df[grep("2-Buten-1-one, 1-\\(2,6,6-trimethyl-1,3-cyclohexadien-1-yl\\)-,.*",cleaned.df$Name),'Name'] <- '2-Buten-1-one, 1-(2,6,6-trimethyl-1,3-cyclohexadien-1-yl)-'

###############################
# Butanoic acid, butyl ester ##
###############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Butanoic acid, butyl ester' & 
    (cleaned.df$X1st.Dimension.Time..s. < 290 | cleaned.df$X1st.Dimension.Time..s. > 291)
),]

###############################
# Butanoic acid, hexyl ester ##
###############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Butanoic acid, hexyl ester' & 
    (cleaned.df$X1st.Dimension.Time..s. < 435 | cleaned.df$X1st.Dimension.Time..s. > 436)
),]


###############################
# à-Farnesene ##
###############################

# removing the (Z,Z)-à-Farnesene
cleaned.df <- cleaned.df[-grep("\\(Z,Z\\)-à-Farnesene",cleaned.df$Name),]

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'à-Farnesene' & 
    (cleaned.df$X1st.Dimension.Time..s. < 670 | cleaned.df$X1st.Dimension.Time..s. > 671)
),]


###################
# Furan, 2-ethyl ##
###################

# fix name
cleaned.df[grep('Furan, 2-ethyl-',cleaned.df$Name),'Name'] <- 'Furan, 2-ethyl'

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Furan, 2-ethyl' & 
    (cleaned.df$X1st.Dimension.Time..s. < 670 | cleaned.df$X1st.Dimension.Time..s. > 671)
),]

###############################
# Butanoic acid, ethyl ester ##
###############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Butanoic acid, ethyl ester' & 
    (cleaned.df$X1st.Dimension.Time..s. < 176 | cleaned.df$X1st.Dimension.Time..s. > 177)
),]

########################################
# Butanoic acid, 2-methylpropyl ester ##
########################################

# fix name
cleaned.df[grep('Butanoic acid, 2-methyl([-][,][ ]){0,}propyl ester',cleaned.df$Name),'Name'] <- 'Butanoic acid, 2-methylpropyl ester'

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Butanoic acid, 2-methylpropyl ester' & 
    (cleaned.df$X1st.Dimension.Time..s. < 236 | cleaned.df$X1st.Dimension.Time..s. > 237)
),]

#####################
# n-Propyl acetate ##
#####################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'n-Propyl acetate' & 
    (cleaned.df$X1st.Dimension.Time..s. < 148 | cleaned.df$X1st.Dimension.Time..s. > 149)
),]

##############################
# Acetic acid, pentyl ester ##
##############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Acetic acid, pentyl ester' & 
    (cleaned.df$X1st.Dimension.Time..s. < 250 | cleaned.df$X1st.Dimension.Time..s. > 300)
),]

############
# Ethanol ##
############

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Ethanol' & 
    (cleaned.df$X1st.Dimension.Time..s. < 131 | cleaned.df$X1st.Dimension.Time..s. > 133)
),]

##############
# Estragole ##
##############

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Estragole' & 
    (cleaned.df$X1st.Dimension.Time..s. < 620 | cleaned.df$X1st.Dimension.Time..s. > 623)
),]

############
# Acetone ##
############

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Acetone' & 
    (cleaned.df$X1st.Dimension.Time..s. < 100 | cleaned.df$X1st.Dimension.Time..s. > 110)
),]

####################
# 4-Oxohex-2-enal ##
####################

first_mean <- mean(cleaned.df[
  which(cleaned.df$Name == '4-Oxohex-2-enal' & 
          (cleaned.df$X1st.Dimension.Time..s. > 500 & cleaned.df$X1st.Dimension.Time..s. < 600)),'X1st.Dimension.Time..s.'])

# fix name for the first compound
cleaned.df[
  which(cleaned.df$Name == '4-Oxohex-2-enal' & 
          (cleaned.df$X1st.Dimension.Time..s. > 500 & 
             cleaned.df$X1st.Dimension.Time..s. < 600)),'Name'] <- '4-Oxohex-2-enal-571'

# fix name for the second compound
cleaned.df[which(cleaned.df$Name == '4-Oxohex-2-enal'),'Name'] <- '4-Oxohex-2-enal-684'

####################
# Tetrahydrofuran ##
####################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Tetrahydrofuran' & 
    (cleaned.df$X1st.Dimension.Time..s. < 116 | cleaned.df$X1st.Dimension.Time..s. > 117)
),]

###############################################
## 1,3,6,10-Dodecatetraene, 3,7,11-trimethyl ##
###############################################

# fix name
cleaned.df[grep("1,3,6,10-Dodecatetraene.*3,7,11-trimethyl.*",cleaned.df$Name),'Name'] <- '1,3,6,10-Dodecatetraene, 3,7,11-trimethyl'

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '1,3,6,10-Dodecatetraene, 3,7,11-trimethyl' & 
    (cleaned.df$X1st.Dimension.Time..s. < 650 | cleaned.df$X1st.Dimension.Time..s. > 660)
),]

##################
# Ethyl Acetate ##
##################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Ethyl Acetate' & 
    (cleaned.df$X1st.Dimension.Time..s. < 119 | cleaned.df$X1st.Dimension.Time..s. >122)
),]

###############
# 2-Heptenal ##
###############

cleaned.df[grep("2-Heptenal,.*",cleaned.df$Name),'Name'] <- '2-Heptenal'

unique(cleaned.df[grep("2-Heptenal.*",cleaned.df$Name),'Name'])

##########################
# 3-Hexen-1-ol, acetate ##
##########################

cleaned.df[grep("3-Hexen-1-ol, acetate,.*",cleaned.df$Name),'Name'] <- '3-Hexen-1-ol, acetate'

########################################
# 5,9-Undecadien-2-one, 6,10-dimethyl ##
########################################

cleaned.df[grep("5,9-Undecadien-2-one, 6,10-dimethyl-.*",cleaned.df$Name),'Name'] <- '5,9-Undecadien-2-one, 6,10-dimethyl'

###################################
# Heptane, 2,2,4,6,6-pentamethyl ##
###################################

# fix name
cleaned.df[which(cleaned.df$Name == 'Heptane, 2,2,4,6,6-pentamethyl-'),'Name'] <- 'Heptane, 2,2,4,6,6-pentamethyl'

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Heptane, 2,2,4,6,6-pentamethyl' & 
    (cleaned.df$X1st.Dimension.Time..s. < 155 | cleaned.df$X1st.Dimension.Time..s. > 157)
),]


##################
# Nitrous oxide ##
##################

# remove nitrous oxide 
cleaned.df <- cleaned.df[-which(cleaned.df$Name == 'Nitrous oxide'),]

##############################
# 2,2,7,7-Tetramethyloctane ##
##############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '2,2,7,7-Tetramethyloctane' & 
    (cleaned.df$X1st.Dimension.Time..s. < 140 | cleaned.df$X1st.Dimension.Time..s. > 150)
),]

####################################
# 1,3-Dioxolane, 2,4,5-trimethyl- ##
####################################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '1,3-Dioxolane, 2,4,5-trimethyl-' & 
    (cleaned.df$X1st.Dimension.Time..s. < 130 | cleaned.df$X1st.Dimension.Time..s. > 150)
),]

#################
# 2-Dodecanone ##
#################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '2-Dodecanone' & 
    (cleaned.df$X1st.Dimension.Time..s. < 640 | cleaned.df$X1st.Dimension.Time..s. > 642)
),]

#############################
# 2(5H)-Furanone, 5-ethyl- ##
#############################

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == '2(5H)-Furanone, 5-ethyl-' & 
    (cleaned.df$X1st.Dimension.Time..s. < 682 | cleaned.df$X1st.Dimension.Time..s. > 685)
),]

#############
# n-Hexane ##
#############

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'n-Hexane' & 
    (cleaned.df$X1st.Dimension.Time..s. < 100 | cleaned.df$X1st.Dimension.Time..s. > 101)
),]

###########
# Nonane ##
###########

# remove noise
cleaned.df <- cleaned.df[-which(
  cleaned.df$Name == 'Nonane' & 
    (cleaned.df$X1st.Dimension.Time..s. < 122 | cleaned.df$X1st.Dimension.Time..s. > 123)
),]

##############
# 2-Butenal ##
##############

# remove the Z configuration because they are only a few (N=27)
cleaned.df <- cleaned.df[-which(cleaned.df$Name == '2-Butenal, (Z)-'),]

# fix the names of other stereoisomers
cleaned.df[which(cleaned.df$Name == '2-Butenal, (E)-'),'Name'] <- '2-Butenal'

write.table(
  cleaned.df,
  'data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated.tsv',
  sep = "\t",
  row.names = FALSE
)

