library(readxl)
library(tidyverse)
gcms_pheno_tbl <- read_excel(
  'data/processed/Supplementary_Data.xlsx',
  sheet = 'GCMS Data'
)

hep_val <- gcms_pheno_tbl$`2-Heptenal`
oct_val <- gcms_pheno_tbl$`(E)-2-Octenal`
fit.hep_oct <- lm(oct_val ~ hep_val)
summary(fit.hep_oct)

data.frame(hep_val,oct_val) %>%
  ggplot(aes(x=hep_val,y=oct_val)) + 
  geom_point(alpha = 0.4, size = 2) +
    geom_line(
      aes(x = hep_val, y = fit.hep_oct$fitted.values),  # predicted data
      color = 'black', alpha = 0.3) +
    GLOBAL_THEME +
    geom_text(aes(label = sprintf("R² = %.2f", summary(fit.hep_oct)$r.squared), x = 2, y = 0.5)) +
    geom_text(aes(label = sprintf("p-value = %s", format(summary(fit.hep_oct)$coefficients[2, 4], scientific = FALSE)), x = 2, y = 0.44)) +
    theme(
      axis.ticks.x = element_blank()
    ) +
    xlab("2-Heptenal") +
    ylab("(E)-2-Octenal")
ggsave(
  filename = 'exploration/2-octenal_and_2-heptenal/2-heptenal_2_octenal_lm.png',
  plot = last_plot(), dpi = 600, width = 5, height=5, units = "in",
  bg = "white"
)
