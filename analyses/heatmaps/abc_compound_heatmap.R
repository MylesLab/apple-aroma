#####################
## LIBRARY IMPORTS ##
#####################

library(tidyverse)
library(readxl)

# load functions
source('themes/theme_main.R')
source('analyses/heatmaps/utils.R')

##################
## DATA LOADING ##
##################

# all the columns that we are interested in for ABC phenotypes
abc_phenotype_names <- c("acidity_17_harv", "percent_acidity_17","brix_17_harv",
"firmness_avg_17_harv", "weight_avg_17_harv", "juiciness_16_harv", "tpc",
"date_jul_17_harv", "flowering_jul_16_harv", "percent_firmness_avg_17")

# get all the ABC phenotypes
abc_pheno_tbl <- 
  read_excel('data/processed/sup_tbl_2-abc_phenotype_table.xlsx') %>% 
  select(all_of(c("apple_id",abc_phenotype_names)))
dim(abc_pheno_tbl)
# [1] 553 11

# get the supplementary table data
gcms_pheno_tbl <- read_excel(
  'data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx'
)
dim(gcms_pheno_tbl)
# [1] 515 107

# load the ABC phenotype data with pretty names
abc_pretty_phenos.df <- read_excel(
  'data/raw/abc_phenotypes_ordered.xlsx'
) %>% select(phenotypes, cleaned)

# load the compound classification data
cmpd_cls.df <- read_excel(
  'data/processed/classification_pivot.xlsx'
)
dim(cmpd_cls.df)
# [1] 106 2

# #########################################################
# ## PLOT ABC PHENOTYPES & COMPOUNDS CORRELATION HEATMAP ##
# #########################################################

# merge the GCMS and ABC phenotype data
pheno_tbl <-
  inner_join(gcms_pheno_tbl, abc_pheno_tbl, by = c("appleid" = "apple_id"))
dim(pheno_tbl)
# [1] 515 117

# create the subset of ABC and GCMS dataframes
pheno_tbl.gcms <- pheno_tbl[,2:107]
pheno_tbl.abc <- pheno_tbl[,108:ncol(pheno_tbl)]

# perform correlation on ABC and GC-MS data
abc_gcms_cor.r <- matrix(NA, nrow = 10, ncol = 106)
abc_gcms_cor.p <- matrix(NA, nrow = 10, ncol = 106)
to_print <- c()
for (abc_pheno_idx in seq_len(ncol(pheno_tbl.abc))) {
  for (gcms_pheno_idx in seq_len(ncol(pheno_tbl.gcms))) {
    abc_pheno  <- as.numeric(unlist(pheno_tbl.abc[, abc_pheno_idx]))
    gcms_pheno <- as.numeric(unlist(pheno_tbl.gcms[, gcms_pheno_idx]))
    test       <- cor.test(abc_pheno, gcms_pheno)
    r          <- as.numeric(test$estimate)
    to_print <- c(
      to_print,
      paste(
        colnames(pheno_tbl.abc)[abc_pheno_idx],
        colnames(pheno_tbl.gcms)[gcms_pheno_idx],
        r
      )
    )
    
    abc_gcms_cor.r[abc_pheno_idx, gcms_pheno_idx] <- r
    abc_gcms_cor.p[abc_pheno_idx, gcms_pheno_idx] <- test$p.value
  }
}

# update the rownames and colnames with the names of variables
rownames(abc_gcms_cor.r) <- colnames(pheno_tbl.abc)
colnames(abc_gcms_cor.r) <- colnames(pheno_tbl.gcms)
rownames(abc_gcms_cor.p) <- colnames(pheno_tbl.abc)
colnames(abc_gcms_cor.p) <- colnames(pheno_tbl.gcms)

# melt the data and create a final dataset with r-value, p-value and the two
# variables
abc_gcms_cor.r.melted <- reshape2::melt(abc_gcms_cor.r) %>% rename(r = value)
abc_gcms_cor.p.melted <- reshape2::melt(abc_gcms_cor.p) %>% rename(p = value)
abc_gcms.cor.melted <- inner_join(abc_gcms_cor.r.melted, abc_gcms_cor.p.melted)
dim(abc_gcms.cor.melted)
# [1] 1060 4

# clean the dataframe
abc_gcms.cor.melted <- 
  left_join(
    abc_gcms.cor.melted, abc_pretty_phenos.df, 
    by = c("Var1" = "phenotypes")
  ) %>%
  rename(abc_phenotype = cleaned, gcms_compound = Var2) %>%
  select(-Var1)

# combine the merged dataset with the class data
abc_gcms.cor.melted <- 
  inner_join(abc_gcms.cor.melted,cmpd_cls.df, by = c("gcms_compound" = "Name"))

abc_gcms.cor.melted %>%
  arrange(Classification, gcms_compound) %>%
  mutate(
    gcms_compound = as_factor(gcms_compound)
  ) %>%
  ggplot(aes(x = gcms_compound, y = abc_phenotype, fill=r)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=10, hjust=1)) +
  coord_fixed() +
  theme(
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    plot.title           = element_text(size = 25, hjust = 0.5, face = "bold"),
    legend.justification = c(1, 0),
    legend.position      = c(0.6, 1),
    legend.direction     = "horizontal"
  ) +
  labs(fill = "Pearson correlation (r)") +
  guides(
    fill = guide_colorbar(
      barwidth       = 10, barheight = 1,
      title.position = "top",
      title.hjust    = 0.5,
      title.theme    = element_text(size = 15, family = "Times")
    )
  )
ggsave(
  filename = "analyses/heatmaps/figures/sup_fig_3a.png",
  plot     = last_plot(),
  dpi      = 600, width = 25, height = 15, units = "in"
)


abc_gcms.cor.melted %>%
  ggplot(aes(x = r)) +
    geom_histogram(colour = "black", fill = "white") +
    scale_y_continuous(expand = c(0, 0)) +
    GLOBAL_THEME +
    xlab("Pearson correlation (r)")
ggsave(
  filename = "analyses/heatmaps/figures/sup_fig_3b.png",
  plot     = last_plot(), bg = "white",
  dpi      = 300, width = 6, height = 6, units = "in"
)

########################
## RESEARCH QUESTIONS ##
########################

bonf_corr_threshold <- 0.05 / nrow(abc_gcms.cor.melted)

## RQ1: How many significant correlations are there?
n_sig_corr <- abc_gcms.cor.melted[abc_gcms.cor.melted$p < bonf_corr_threshold,] %>% nrow()
n_sig_corr
# [1] 69
nrow(abc_gcms.cor.melted)
# [1] 1060
# Of the total of 1,060 comparisons, 69 (~6.5%) are significant.

## RQ2: How many +'ve significant correlations are there?
n_pos_sig_corr <- abc_gcms.cor.melted[
  (abc_gcms.cor.melted$p < bonf_corr_threshold) & (abc_gcms.cor.melted$r > 0),
] %>% nrow()
n_pos_sig_corr
# [1] 9
# There are 9 (~0.8% of all comparisons and ~13% of all the significant correlations) +'ve significant correlations.

## RQ2: How many -'ve significant correlations are there?
n_neg_sig_corr <- abc_gcms.cor.melted[
  (abc_gcms.cor.melted$p < bonf_corr_threshold) & (abc_gcms.cor.melted$r < 0),
] %>% nrow()
n_neg_sig_corr
# [1] 60
# There are 60 (~5.7% of all comparisions and ~87% of all significant correlations) -'ve significant correaltions

## RQ3: Are there more positive correlations than there are negative?
# we will answer this through a chi-squared test
chisq.test(as.table(rbind(
  Observed = c(
    Positive = n_pos_sig_corr,
    Negative = n_neg_sig_corr
  ),
  Expected = c(
    Positive = (n_sig_corr / 2),
    Negative = (n_sig_corr / 2)
  )
)))
# p-value: 7.157396e-06
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  as.table(rbind(Observed = c(Positive = n_pos_sig_corr, Negative = n_neg_sig_corr),     Expected = c(Positive = (n_sig_corr/2), Negative = (n_sig_corr/2))))
# X-squared = 20.151, df = 1, p-value = 7.157e-06



