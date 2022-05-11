# CONTEXT: Ardra Bio is interested in exploring the potential of trans-2-hexenal (green leaf volatile compound) for the
#          preservation of fruits and vegetables. They are hoping to apply for AAFC grant for the category of
#          "technologies that can help reduce food waste".
#
# We need to investigate the correlation of our phenotypes with 2-hexanal, to see if there are any strong correlations.
# These correlations will indicate an importance of this compound in the particular phenotypes.

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

source('themes/theme_main.R')

generate_dot_plot <- function(x, y, y_name, r2_val, pval, positions) {

  x1 <- as.numeric(positions[1])
  x2 <- as.numeric(positions[2])
  y1 <- as.numeric(positions[3])
  y2 <- as.numeric(positions[4])

  ggplot(figures.df, aes_string(x, y)) +
    stat_summary(fun.data = mean_cl_normal, alpha = 0.5) +
    geom_smooth(method = 'lm', show.legend = T, formula = y ~ x) +
    geom_text(x = x1, y = y1, label = r2_val) +
    geom_text(x = x2, y = y2, label = pval) +
    GLOBAL_THEME +
    xlab("Hexanal") +
    ylab(y_name)
}

# load the GC-MS dataset
gcms_data.tbl <- read_excel('data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx')
dim(gcms_data.tbl)
# [1] 515 107

# gcms_data.noaid.tbl <- gcms_data.tbl[,2:ncol(gcms_data.tbl)]

# load the phenotype dataset
abc_pheno.tbl <- read_excel('data/processed/sup_tbl_2-abc_phenotype_table_v2.xlsx')
dim(abc_pheno.tbl)
# [1] 515 44

# merge the two data togther
final.df <- inner_join(gcms_data.tbl, abc_pheno.tbl, by = c("appleid" = "apple_id"))
dim(final.df)
# [1] 515 150


summary(final.df$`2-Hexenal`)

# generating a dataframe for figures
figures.df <- data.frame(
  TwoHexenal = final.df$Hexanal,
  DeltaFirmness = final.df$percent_firmness_avg_17,
  DeltaAcidity = final.df$percent_acidity_17,
  Weight = final.df$weight_avg_17_harv,
  PhenolicContent = final.df$tpc,
  HarvestDate = final.df$date_jul_17_harv,
  FloweringTime = final.df$flowering_jul_16_harv,
  Firmness = final.df$firmness_avg_17_harv,
  Brix = final.df$brix_17_harv,
  Acidity = final.df$acidity_17_harv
)

# list of phenotype and their relevant data.
# here's the description of the columns:
# 1. the name of the column in df
# 2. Pretty name of the phenotype
# 3. x position for r2 value
# 4. x position for pvalue
# 5. y position for r2 value
# 6. y position for pvalue
phenotypes <- list(
  c("DeltaFirmness", "Δ Firmness (%)", 60, 60, -9, -14),
  c("DeltaAcidity", "Δ Acidity (%)", 60, 60, 2, -4),
  c("Weight", "Weight", 60, 60, 312, 287),
  c("PhenolicContent", "PhenolicContent", 40, 40, 22, 19),
  c("HarvestDate", "Harvest Date (Julian days)", 60, 60, 240, 236),
  c("FloweringTime", "Flowering Time (Julian days)", 60, 60, 152, 149),
  c("Firmness", "Firmness (%)", 60, 60, 6.5, 5.5),
  c("Brix", "Δ Brix (%)", 60, 60, 15.5, 14.5),
  c("Acidity", "Acidity", 60, 60, 22, 19)
)

plots <- list()
for (p in seq_along(phenotypes)) {
  pheno <- phenotypes[[p]][1]
  pretty_name <- phenotypes[[p]][2]

  x1 <- phenotypes[[p]][3]
  x2 <- phenotypes[[p]][4]
  y1 <- phenotypes[[p]][5]
  y2 <- phenotypes[[p]][6]

  m <- lm(get(pheno, figures.df) ~ figures.df$TwoHexenal)
  r2 <- round(summary(m)$r.squared, 3)
  pval <- format(summary(m)$coefficients[2, 4], scientific = TRUE)

  plots[[p]] <- generate_dot_plot(
    "TwoHexenal",
    phenotypes[[p]][1],
    pretty_name,
    paste0("R² = ", r2),
    paste0("P-value = ", pval),
    c(x1, x2, y1, y2)
  )

}
figures.png <- ggarrange(plotlist = plots, ncol = 3, nrow = 3)
ggsave(
  'exploration/explore-correlation-of-2-hexenal-with-abc-phenotypes/correlation_of_2-hexenal_with_abc_phenotypes.png',
  figures.png,
  width = 1412,
  height = 640,
  dpi = 80,
  units = "px",
  bg = "white"
)