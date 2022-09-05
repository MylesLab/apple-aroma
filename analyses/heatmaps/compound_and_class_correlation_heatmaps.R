#####################
## LIBRARY IMPORTS ##
#####################

library(tidyverse)
library(readxl)
library(ggpubr)

# load functions
source('themes/theme_main.R')
source('analyses/heatmaps/utils.R')

##################
## DATA LOADING ##
##################

# get the supplementary table data
gcms_pheno_tbl <- read_excel(
  'data/processed/Supplementary_Data.xlsx',
  sheet = 'GCMS Data'
)
dim(gcms_pheno_tbl)
# [1] 515 107

gcms_pheno_tbl.noaid <- gcms_pheno_tbl[, 2:ncol(gcms_pheno_tbl)]
dim(gcms_pheno_tbl.noaid)
# [1] 515 106

# load the compound classification data
cmpd_cls.df <- read_excel(
  'data/processed/Supplementary_Data.xlsx',
  sheet = 'Compound Class'
)
dim(cmpd_cls.df)
# [1] 106 2

##################################
## PERFORM CORRELATION ANALYSIS ##
##################################

# perform the correlation
gcms_data.cor.df <- perform_matrix_correlation(gcms_pheno_tbl.noaid)

gcms_data.cor.r <- gcms_data.cor.df$r
dim(gcms_data.cor.r)
# [1] 106 106 

gcms_data.cor.pval <- gcms_data.cor.df$pval
dim(gcms_data.cor.pval)
# [1] 106 106

# melt the correlation r and p matrices
gcms_data.cor.r.melted <- reshape2::melt(gcms_data.cor.r) %>%
  rename(r = value)
dim(gcms_data.cor.r.melted)
# [1] 11236 3

gcms_data.cor.p.melted <- reshape2::melt(gcms_data.cor.pval) %>%
  rename(p = value)
dim(gcms_data.cor.p.melted)
# [1] 11236 3

# combine the r and p into one data frame
gcms_data.cor.melted <- inner_join(gcms_data.cor.r.melted, gcms_data.cor.p.melted)
dim(gcms_data.cor.melted)
# [1] 11236 4

# combine the compound classification data with the correlation data
gcms_data.cor.melted.cls <- left_join(gcms_data.cor.melted, cmpd_cls.df, by = c("Var1" = "Compound name"))
dim(gcms_data.cor.melted.cls)
# [1] 11236 5

gcms_data.cor.melted.cls <- left_join(gcms_data.cor.melted.cls, cmpd_cls.df, by = c("Var2" = "Compound name"))
dim(gcms_data.cor.melted.cls)
# [1] 11236     7

#######################################
## PLOT COMPOUND CORRELATION HEATMAP ##
#######################################

## PLOTTING COMPOUND CORRELATION HEATMAP

# create a dataframe for plotting Figure 2
sup_fig_2.df <- gcms_data.cor.melted.cls %>%
  arrange(Classification.x) %>%
  mutate(
    Var1 = as_factor(Var1),
    Var2 = factor(Var2, levels = levels(Var1))
  ) %>%
  filter(as.integer(Var1) > as.integer(Var2))

dim(sup_fig_2.df)
# [1] 5565 6

# save the figure 2A data
openxlsx::write.xlsx(
  sup_fig_2.df,
  file  = "data/processed/figures/sup_fig_2.xlsx",
)

# plot the full-length correlation heatmap
sup_fig_2a.plot <- sup_fig_2.df %>%
  ggplot(aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1, 1)) +
  theme_classic() +
  theme(
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    axis.text.x          = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.6, 0.7),
    legend.direction     = "horizontal",
    plot.title           = element_text(size = 25, hjust = 0.5, face = "bold")
  ) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE), position = "right") +
  coord_fixed() +
  labs(fill = "Pearson correlation (r)") +
  guides(
    fill = guide_colorbar(
      barwidth       = 30, barheight = 2,
      title.position = "top",
      title.hjust    = 0.5,
      title.theme    = element_text(size = 34, family = "Times")
    )
  )

# save the plot
ggsave(
  filename = "analyses/heatmaps/figures/sup_fig_2a.png",
  plot     = sup_fig_2a.plot,
  bg       = "white",
  dpi      = 600, width = 25, height = 15, units = "in"
)

## PLOTTING DISTRIBUTION OF r-VALUES
sup_fig_2b.plot <- sup_fig_2.df %>%
  ggplot(aes(x = r)) +
  geom_histogram(colour = "black", fill = "white") +
  scale_y_continuous(expand = c(0, 0)) +
  GLOBAL_THEME +
  xlab("Pearson correlation (r)")

ggsave(
  filename = "analyses/heatmaps/figures/sup_fig_2b.png",
  plot     = sup_fig_2b.plot, bg = "white",
  dpi      = 300, width = 6, height = 6, units = "in"
)

bonf_corr_thresh <- 0.05 / ((106 * 105) / 2)

## PLOTTING DISTRIBUTION OF p-VALUES
sup_fig_2c.plot <- sup_fig_2.df %>%
  ggplot(aes(x = -log10(p))) +
  geom_histogram(colour = "black", fill = "white", bins = 40) +
  geom_vline(xintercept = bonf_corr_thresh, colour = "red") +
  scale_y_continuous(expand = c(0, 0)) +
  GLOBAL_THEME +
  xlab("-log10(p)")

ggsave(
  filename = "analyses/heatmaps/figures/sup_fig_2c.png",
  plot     = sup_fig_2c.plot, bg = "white",
  dpi      = 300, width = 6, height = 6, units = "in"
)

## RQ1: How many significant correlations are there?
n_sig_corr <- sup_fig_2.df[sup_fig_2.df$p < bonf_corr_thresh,] %>% nrow()
# [1] 726
nrow(sup_fig_2.df)
# [1] 5565
# Of the total of 5,565 comparisons, 726 (~13%) are significant.

## RQ2: How many +'ve significant correlations are there?
pos_sig_corr.df <- sup_fig_2.df[
  (sup_fig_2.df$p < bonf_corr_thresh) & (sup_fig_2.df$r > 0),
]
n_pos_sig_corr <- nrow(pos_sig_corr.df)
# [1] 704
# There are 704 (~13% of all comparisons and ~96% of all the significant correlations) +'ve significant correlations.

# save the top 50 positive correlations into a file
pos_sig_corr.df %>%
  rename(
    "Classification Var1" = "Classification.x",
    "Classification Var2" = "Classification.y"
  ) %>% arrange(desc(r)) %>%
  select("Var1", "Var2", "Classification Var1", "Classification Var2") %>%
  head(50) %>% 
  openxlsx::write.xlsx(
    .,
    file = "data/processed/heatmaps/top_50_positive_significant_correlations_by_r-value.xlsx")


## RQ2: How many -'ve significant correlations are there?
neg_sig_corr.df <- sup_fig_2.df[
  (sup_fig_2.df$p < bonf_corr_thresh) & (sup_fig_2.df$r < 0),
]
n_neg_sig_corr <- nrow(neg_sig_corr.df)
# [1] 22
# There are 22 (~0.4% of all comparisions and ~3% of all significant correlations) -'ve significant correaltions

# save the negative correlations into a file
neg_sig_corr.df %>%
  rename(
    "Classification Var1" = "Classification.x",
    "Classification Var2" = "Classification.y") %>%
  select("Var1", "Var2", "Classification Var1", "Classification Var2") %>%
  openxlsx::write.xlsx(
    .,
    file = "data/processed/heatmaps/negative_significant_correlations.xlsx"
)

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
# p-value: 7.390595e-91
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  as.table(rbind(Observed = c(Positive = n_pos_sig_corr, Negative = n_neg_sig_corr),     Expected = c(Positive = (n_sig_corr/2), Negative = (n_sig_corr/2))))
# X-squared = 408.6, df = 1, p-value < 2.2e-16

# There are statistically significantly more positive correlations than there 
# are negative correlations.

## RQ4: Are significant +'ve correlations more often within class than between class?
n_within_class_comp <- sup_fig_2.df[
  (sup_fig_2.df$Classification.x == sup_fig_2.df$Classification.y),
] %>% nrow()
n_within_class_comp
# [1] 764

n_between_class_comp <- sup_fig_2.df[
  (sup_fig_2.df$Classification.x != sup_fig_2.df$Classification.y),
] %>% nrow()
n_between_class_comp
# [1] 4801


n_within_class_sig_comp <- sup_fig_2.df[
  (sup_fig_2.df$p < bonf_corr_thresh) & (sup_fig_2.df$Classification.x == sup_fig_2.df$Classification.y),
] %>% nrow()
n_within_class_sig_comp
# [1] 235
# There are 235 (4.2% of all comparisons; 32% of sig. comparisons) that are within the same class

n_between_class_sig_comp <- sup_fig_2.df[
  (sup_fig_2.df$p < bonf_corr_thresh) & (sup_fig_2.df$Classification.x != sup_fig_2.df$Classification.y),
] %>% nrow()
n_between_class_sig_comp
# [1] 491
# There are 491 (9% of all comparisons; 68% of sig. comparisons) that are between different classes

chisq.test(as.table(rbind(
  Observed = c(
    Within = n_within_class_sig_comp,
    Between = n_between_class_sig_comp
  ),
  Expected = c(
    Within = (n_within_class_comp / nrow(sup_fig_2.df) ) * n_sig_corr,
    Between = (n_between_class_comp / nrow(sup_fig_2.df) ) * n_sig_corr
  )
)))
# p-value: 5.732265e-17
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  as.table(rbind(Observed = c(Within = n_within_class_sig_comp,     Between = n_between_class_sig_comp), Expected = c(Within = (n_within_class_comp/5565) *     n_sig_corr, Between = (n_between_class_comp/5565) * n_sig_corr)))
# X-squared = 70.067, df = 1, p-value < 2.2e-16


# ######################################
# ## PLOT CLASSES CORRELATION HEATMAP ##
# ######################################

# create a data frame with cumulative sum of all the compounds based on the classes.
new_colnames                   <- cmpd_cls.df[match(colnames(gcms_pheno_tbl.noaid), cmpd_cls.df$`Compound name`),]$Classification
gcms_class_tbl.noaid           <- gcms_pheno_tbl.noaid
colnames(gcms_class_tbl.noaid) <- new_colnames

# sum up all the abundance value for each class for each sample
gcms_cls.df <- rowsum(
  t(gcms_class_tbl.noaid),
  group = colnames(gcms_class_tbl.noaid),
  na.rm = TRUE
) %>% t()

gcms_cls_cor.df <- perform_matrix_correlation(gcms_cls.df)

gcms_cls.cor.r <- gcms_cls_cor.df$r
gcms_cls.cor.pval <- gcms_cls_cor.df$pval

# melt the correlation r and p matrices
gcms_cls.cor.r.melted <- reshape2::melt(gcms_cls.cor.r) %>%
  rename(r = value)
dim(gcms_cls.cor.r.melted)
# [1] 169 3

gcms_cls.cor.p.melted <- reshape2::melt(gcms_cls.cor.pval) %>%
  rename(p = value)
dim(gcms_cls.cor.p.melted)
# [1] 169 3


# combine the r and p into one data frame
gcms_cls.cor.melted <- inner_join(gcms_cls.cor.r.melted, gcms_cls.cor.p.melted)
dim(gcms_cls.cor.melted)
# [1] 169 4

# only keep the lower triangle values
gcms_cls.cor.melted <- gcms_cls.cor.melted %>%
  filter(as.integer(Var1) > as.integer(Var2))
dim(gcms_cls.cor.melted)
# [1] 78 4

# plot the full-length correlation heatmap
fig_2a.plot <- gcms_cls.cor.melted %>%
  ggplot(aes(x = Var1, y = Var2, fill = r)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1, 1)) +
  theme_classic() +
  theme(
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    axis.text.x          = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    panel.background     = element_blank(),
    axis.ticks           = element_blank(),
    legend.justification = c(1, 0),
    legend.position      = c(0.6, 0.7),
    legend.direction     = "horizontal",
    plot.title           = element_text(size = 15, hjust = 0.5, face = "bold")
  ) +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  scale_y_discrete(guide = guide_axis(check.overlap = TRUE), position = "right") +
  coord_fixed() +
  labs(fill = "Pearson correlation (r)") +
  guides(
    fill = guide_colorbar(
      barwidth       = 15, barheight = 1,
      title.position = "top",
      title.hjust    = 0.5,
      title.theme    = element_text(size = 14, family = "Times")
    )
  )

ggsave(
  filename = "analyses/heatmaps/figures/fig_2a.png",
  plot     = fig_2a.plot,
  bg       = "white",
  dpi      = 300, width = 15, height = 5, units = "in"
)

## PLOTTING DISTRIBUTION OF r-VALUES
fig_2b.plot <- gcms_cls.cor.melted %>%
  ggplot(aes(x = r)) +
  geom_histogram(colour = "black", fill = "white") +
  scale_y_continuous(expand = c(0, 0)) +
  GLOBAL_THEME +
  xlab("Pearson correlation (r)")

ggsave(
  filename = "analyses/heatmaps/figures/fig_2b.png",
  plot     = fig_2b.plot, bg = "white",
  dpi      = 300, width = 6, height = 6, units = "in"
)

fig_2_bonf_corr_thresh <- 0.05 / ((13 * 12) / 2)

## PLOTTING DISTRIBUTION OF p-VALUES
fig_2c.plot <- gcms_cls.cor.melted %>%
  ggplot(aes(x = -log10(p))) +
  geom_histogram(colour = "black", fill = "white", bins = 40) +
  geom_vline(xintercept = -log10(fig_2_bonf_corr_thresh), colour = "red") +
  scale_y_continuous(expand = c(0, 0)) +
  GLOBAL_THEME +
  xlab("-log10(p)")

ggsave(
  filename = "analyses/heatmaps/figures/fig_2c.png",
  plot     = fig_2c.plot, bg = "white",
  dpi      = 300, width = 6, height = 6, units = "in"
)


fig2_plots      <- list()
fig2_plots[[1]] <- fig_2a.plot
fig2_plots[[2]] <- ggarrange(
  fig_2b.plot, fig_2c.plot,
  nrow          = 1, ncol = 2,
  common.legend = T,
  labels        = c("B", "C")
)
figure_2.plot   <- ggarrange(plotlist = fig2_plots, nrow = 1, ncol = 2, labels = c("A", ""))
ggsave(
  filename = "figures/final_figures/Figure_2.png",
  figure_2.plot,
  dpi      = 600,
  bg       = "white",
  width    = 20,
  height   = 5,
  units    = "in"
)
