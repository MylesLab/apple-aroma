library(dplyr)

source('exploration/genomic-inflation-calculation/utils.R')
source('themes/theme_main.R')

#######################################################################
## GENERATE THE LAMBDA VECTOR AND SAVE MANHATTAN AND QQ-PLOT FIGURES ##
#######################################################################

## FOR RAW GWAS
tictoc::tic("RAW GWAS")
raw_gwas_lambdas.df <- generate_gi_lambda_data(
  "/Users/tayabsoomro/Documents/Projects/GCMS-Project/gwas_1",
  "exploration/genomic-inflation-calculation/figures_gwas_1/",
  alert_when_done = TRUE
)
tictoc::toc()
beepr::beep(sound=3)
dim(raw_gwas_lambdas.df)
# [1] 104 2

## FOR CLASS-CORRECTED GWAS
tictoc::tic("CLASS CORRECTED GWAS")
class_corrected_gwas_lambdas.df <- generate_gi_lambda_data(
  "/Users/tayabsoomro/Documents/Projects/GCMS-Project/class_corrected_gwas_analysis/gwas_data_updated_names/all-samples",
  "exploration/genomic-inflation-calculation/figures_class_corrected_gwas/",
  alert_when_done = TRUE
)
tictoc::toc()
beepr::beep(sound=3)

tictoc::tic("CLASS GWAS")
class_gwas_lambdas.df <- generate_gi_lambda_data(
  "/Users/tayabsoomro/Documents/Projects/GCMS-Project/class_gwas_analysis/gwas_data/all-samples",
  "exploration/genomic-inflation-calculation/figures_class_gwas/",
  alert_when_done = TRUE
)
tictoc::toc()

############################################
## GENERATE THE SHAPIRO-WILK TEST RESULTS ##
############################################

raw_gwas_wilk.df <- generate_wilcoxon_data(
  "data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx"
)

##########################
## CREATE FINAL DATASET ##
##########################

final.df <- left_join(
  raw_gwas_lambdas.df,
  raw_gwas_wilk.df
)

final.df$logp <- -log10(as.numeric(final.df$p))


###############
## VISUALIZE ##
###############

# W-statistic
final.df %>% ggplot(aes(x = W, y = Lambda)) +
  geom_point(alpha = 0.4, size = 2, color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "green") +
  GLOBAL_THEME +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  ylab(expression("Genomic Inflation Factor (  " * lambda * " )")) +
  xlab("W-statistic")

# p-value
final.df %>% ggplot(aes(x = logp, y = Lambda)) +
  geom_point(alpha = 0.4, size = 2, color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "green") +
  geom_vline(xintercept = -log10(0.05 / nrow(final.df)), linetype = "dashed", color = "red") +
  GLOBAL_THEME +
  theme(
    axis.ticks.x = element_blank()
  ) +
  ylab(expression("Genomic Inflation Factor (  " * lambda * " )")) +
  xlab("-log10(p)")


############################################
## GENERATE THE SAMPLE MISSINGNESS VECTOR ##
############################################

# load the genotype table
gcms_pheno_tbl <- read_excel(
  'data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx'
)
dim(gcms_pheno_tbl)
# [1] 515 107

# removing the two compounds for for which the GWAS didn't run
subset_gcms_pheno_tbl <- gcms_pheno_tbl[, !names(gcms_pheno_tbl) %in% c("1-Pentanol", "Pentanal", "appleid")]
dim(subset_gcms_pheno_tbl)
# [1] 515 104

# create a data frame of missingness and compound
gcms_missingness_df <- as.data.frame(colSums(subset_gcms_pheno_tbl == 0))
gcms_missingness_df$compound <- rownames(gcms_missingness_df)
rownames(gcms_missingness_df) <- NULL
colnames(gcms_missingness_df) <- c("num_missing", "compound")
gcms_missingness_df$percent_missing <- (gcms_missingness_df$num_missing / 515) * 100

# combine the missingess and lambda columns
final.df <- left_join(lambdas, gcms_missingness_df, by = "compound")

fit_lambda_missing <- lm(lambda ~ percent_missing, final.df)
summary(fit_lambda_missing)

cor.test(final.df$lambda, final.df$percent_missing, method = "spearman")
summary(fit_lambda_missing)
final.df %>% ggplot(aes(x = percent_missing, y = lambda)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_line(
    aes(x = percent_missing, y = fit_lambda_missing$fitted.values),  # predicted data
    color = 'black', alpha = 0.3) +
  geom_hline(yintercept = 1, color = "red") +
  geom_text(aes(label = sprintf("R² = %f", summary(fit_lambda_missing)$r.squared), x = 25, y = 0.7)) +
  geom_text(aes(label = sprintf("p-value = %s", format(summary(fit_lambda_missing)$coefficients[2, 4], scientific = FALSE)), x = 25, y = 0.66)) +
  GLOBAL_THEME +
  theme(
    axis.ticks.x = element_blank()
  ) +
  xlab("Percent Missing Compounds (%)") +
  ylab(expression("Genomic Inflation Factor (" * lambda * ")"))


View(final.df[which(final.df$lambda > 1.5),])