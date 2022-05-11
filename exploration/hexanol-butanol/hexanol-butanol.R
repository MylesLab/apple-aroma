library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

source('themes/theme_main.R')

# load the GCMS table
gcms_pheno_tbl <- read_excel(
  'data/processed/Supplementary_Data.xlsx',
  sheet = 'GCMS Data'
)
dim(gcms_pheno_tbl)
# [1] 515 107

# get the volatile data for both compounds
butanol_dat <- gcms_pheno_tbl[, c("appleid", "1-Butanol")]
hexanol_dat <- gcms_pheno_tbl[, c("appleid", "1-Hexanol")]

# load the geno type data
genotype_data <- read.table(
  'data/raw/nac_locus.tsv',
  sep = " ",
  header = TRUE
)


# get individual apple ids for different genotypes
ind_0 <- data.frame(
  appleid = genotype_data[genotype_data$X3_30698039_vineland_nac_A == 0, 'FID'],
  genotype = "CC"
)
ind_1 <- data.frame(
  appleid = genotype_data[genotype_data$X3_30698039_vineland_nac_A == 1, 'FID'],
  genotype = "AC"
)
ind_2 <- data.frame(
  appleid = genotype_data[genotype_data$X3_30698039_vineland_nac_A == 2, 'FID'],
  genotype = "AA"
)

# get the 1-Butanol volatile data for each genotype
butanol_0 <- left_join(ind_0, butanol_dat)
butanol_1 <- left_join(ind_1, butanol_dat)
butanol_2 <- left_join(ind_2, butanol_dat)
final.butanol.df <- rbind(butanol_0, butanol_1, butanol_2)

# get the 1-Hexanol volatile data for each genotype
hexanol_0 <- left_join(ind_0, hexanol_dat)
hexanol_1 <- left_join(ind_1, hexanol_dat)
hexanol_2 <- left_join(ind_2, hexanol_dat)
final.hexanol.df <- rbind(hexanol_0, hexanol_1, hexanol_2)

# combine everything for plotting
final.df <- left_join(final.butanol.df, final.hexanol.df)

melted_final.df <- reshape2::melt(final.df, id.vars = c("appleid", "genotype"), measure.vars = c("1-Butanol", "1-Hexanol"))

num_AA <- melted_final.df %>%
  filter(genotype == "AA") %>%
  nrow()
num_AC <- melted_final.df %>%
  filter(genotype == "AC") %>%
  nrow()
num_CC <- melted_final.df %>%
  filter(genotype == "CC") %>%
  nrow()

GLOBAL_LABELS <- c(
  paste0("N=", num_AA),
  paste0("N=", num_AC),
  paste0("N=", num_CC)
)

nac_region_genotypes.plot <-
  ggplot(melted_final.df, aes(x = genotype, y = value, fill = variable)) +
    geom_boxplot() +
    GLOBAL_THEME +
    xlab("") +
    ylab("Total Ion \nCount (TIC)") +
    labs(fill = "") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 28),
      axis.text.y = element_text(size = 28),
      axis.title.x = element_text(size = 28),
      axis.title.y = element_text(size = 28),
      legend.text = element_text(size = 28)
    ) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    annotate(geom = 'text', y = -5, x = 1, label = paste0("N=", num_AA)) +
    annotate(geom = 'text', y = -5, x = 2, label = paste0("N=", num_AC)) +
    annotate(geom = 'text', y = -5, x = 3, label = paste0("N=", num_CC))

ggsave(
  filename = "exploration/hexanol-butanol/figures/nac_region_hexanol_butanol_genotypes.png",
  nac_region_genotypes.plot,
  width = 815,
  height = 515,
  dpi = 75,
  units = "px"
)
nac_region_genotypes.plot

# linear relationship between 1-Butanol and 1-Hexanol
fit.butanol_vs_hexanol <- lm(final.df$`1-Butanol` ~ final.df$`1-Hexanol`)
summary(fit.butanol_vs_hexanol)
ggplot(final.df, aes(x = `1-Hexanol`, y = `1-Butanol`)) +
  geom_point(alpha = 0.1, size = 7) +
  geom_line(
    aes(x = `1-Hexanol`, y = fit.butanol_vs_hexanol$fitted.values),  # predicted data
    color = 'black', alpha = 0.3) +
  geom_text(size = 5, aes(label = sprintf("R² = %f", summary(fit.butanol_vs_hexanol)$r.squared), x = 85, y = 10)) +
  geom_text(size = 5, aes(label = sprintf("p-value = %s", format(summary(fit.butanol_vs_hexanol)$coefficients[2, 4], scientific = TRUE)), x = 85, y = 8)) +
  GLOBAL_THEME +
  theme(
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )
ggsave(
  filename = "exploration/hexanol-butanol/figures/hexanol-butanol_regression.png",
  plot = last_plot(),
  bg = "white"
)
