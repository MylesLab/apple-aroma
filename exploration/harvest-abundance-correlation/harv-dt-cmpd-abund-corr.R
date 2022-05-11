############################
## IMPORTS & DATA LOADING ##
############################

library(readxl)
library(dplyr)

# load the theme for plots
source("themes/theme_main.R")

###############################
## ORGANIZE THE LARSEN DATA ##
##############################

larsen_gcms_dat <- as.data.frame(read_excel('data/raw/larsen/larsen_data.xlsx'))
dim(larsen_gcms_dat)
# [1] 145 51
# There seems to be 145 cultivars and 49 compounds detected, this is consistent with the paper

# calculating the abundance of these compounds
compound_col_idx <- 3:ncol(larsen_gcms_dat)
larsen_gcms_dat$CumulativeAbundance <- rowSums(larsen_gcms_dat[, compound_col_idx])

# load the harvest date data
larsen_harv_dat <- as.data.frame(read_excel('data/raw/larsen/larsen_data_harvest_date_info.xlsx'))
dim(larsen_harv_dat)
# [1] 177 3

# remove the rows for which harvest date is "NA"
larsen_harv_dat <- larsen_harv_dat[!(larsen_harv_dat$`Harvest date` %in% "NA"),]

dim(larsen_harv_dat)
# [1] 159 3

# combine the harvest date and GCMS datasets together
larsen_combined_dat <- inner_join(larsen_gcms_dat, larsen_harv_dat, by = "Acc. Nr.")
dim(larsen_combined_dat)
# [1] 122 54

# We retain 122 cultivars which are present in both the harvest date and
# GCMS data.

# checking to see if the joining was successful
sum(larsen_harv_dat$`Acc. Nr.` %in% larsen_gcms_dat$`Acc. Nr.`)
# [1] 122

sum(larsen_harv_dat$`Accession name` %in% larsen_gcms_dat$`Accession name`)
# [1] 122

#######################
## ORGANIZE OUR DATA ##
#######################

# loading the ABC phenotype table
abc_pheno_tbl <- read_excel('data/processed/sup_tbl_2-abc_phenotype_table_v2.xlsx')
dim(abc_pheno_tbl)
# [1] 515 44

# loading the GCMS phenotype table
gcms_pheno_tbl <- read_excel(
  'data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx')
dim(gcms_pheno_tbl)
# [1] 515 107

# join the ABC and GCMS data so that we can obtain the harvest date and compound abundance
our_combined_dat <- inner_join(gcms_pheno_tbl, abc_pheno_tbl, by = c("appleid" = "apple_id"))
dim(our_combined_dat)
# [1] 515 150

# check to see if the joining was successful
sum(abc_pheno_tbl$apple_id %in% gcms_pheno_tbl$appleid)
# [1] 515
# great, the numbers match

aroma_columns_idx <- seq_along(gcms_pheno_tbl)[-1]

# calculate the cumulative sum of all the compound concentrations
our_combined_dat$CumulativeAbundance <- rowSums(our_combined_dat[, aroma_columns_idx])

###############
## FUNCTIONS ##
###############

hrv_dt_abund_corr_lm <- function(harvest_date_dat, cumulative_abund_dat) {
  fit.hv_missing <- lm(cumulative_abund_dat ~ harvest_date_dat)
  summary(fit.hv_missing)
  plt <- ggplot(data.frame(cumulative_abund_dat, harvest_date_dat), aes(x = harvest_date_dat, y = cumulative_abund_dat)) +
    geom_point(alpha = 0.5, color = "black") +                           # observed data
    geom_line(aes(x = harvest_date_dat, y = fit.hv_missing$fitted.values),  # predicted data
              color = 'black', alpha = 0.3) +
    GLOBAL_THEME +
    xlab("Harvest Date") +
    ylab("Cumulative Compound Abundance")

  return(list(ModelResult = fit.hv_missing, Plot = plt))
}

########################
## LARSEN CORRELATION ##
########################

larsen_result <- hrv_dt_abund_corr_lm(
  as.numeric(larsen_combined_dat$`Harvest date`),
  as.numeric(larsen_combined_dat$CumulativeAbundance)
)
summary(larsen_result$ModelResult)
# Call:
#   lm(formula = cumulative_abund_dat ~ harvest_date_dat)
#
# Residuals:
#   Min         1Q     Median         3Q        Max
# -297436452 -141240545  -30183729  113548779  522301265
#
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)      328782914   45734887   7.189 6.05e-11 ***
# harvest_date_dat  -8495840   10530990  -0.807    0.421
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 175100000 on 120 degrees of freedom
# Multiple R-squared:  0.005394,	Adjusted R-squared:  -0.002894
# F-statistic: 0.6508 on 1 and 120 DF,  p-value: 0.4214

larsen_result$Plot

##########################
## OUR GCMS CORRELATION ##
##########################

our_result <- hrv_dt_abund_corr_lm(
  our_combined_dat$date_jul_17_harv,
  our_combined_dat$CumulativeAbundance
)
summary(our_result$ModelResult)
# Call:
#   lm(formula = cumulative_abund_dat ~ harvest_date_dat)
#
# Residuals:
#   Min      1Q  Median      3Q     Max
# -218.34 -106.67  -36.25   85.55  625.56
#
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)       719.580     96.556   7.452 3.92e-13 ***
# harvest_date_dat   -1.918      0.367  -5.226 2.52e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 136 on 513 degrees of freedom
# Multiple R-squared:  0.05055,	Adjusted R-squared:  0.0487
# F-statistic: 27.31 on 1 and 513 DF,  p-value: 2.524e-07

our_result$Plot

############
# RESULTS ##
############