################
## FIGS INDEX ##
################

# FIGS
## fig_1a - Distribution of total abundance & ubiquity
## fig_1d - Distribution of volatiles detected and abundance
## fig_1e - Distribution of volatiles abundance
## fig_1b - Distribution of volatile classes detected
## fig_1c - Distribution of volatile classes abundance
## fig_s1 - Distribution of Sample Ubiquity for Volatiles


#####################
## LIBRARY IMPORTS ##
#####################

library(readxl)
library(tidyverse)
library(ggpubr)
source("themes/theme_main.R")
source("utils/basic_stats.R")

library(grid)
library(forcats)
library("viridis")

#############################
## DATA LOADING & CURATION ##
#############################

# loading the GCMS phenotype table
gcms_pheno_tbl <- read_excel(
  "data/processed/Supplementary_Data.xlsx",
  sheet = "GCMS Data"
)
dim(gcms_pheno_tbl)
# [1] 515 107

# create a GCMS phenotype table without apple id
gcms_pheno_noaid_tbl <- gcms_pheno_tbl[, 2:ncol(gcms_pheno_tbl)]
dim(gcms_pheno_noaid_tbl)
# [1] 515 106

# load the volatile classification data
classification_pivot_tbl <- read_excel(
  "data/processed/Supplementary_Data.xlsx",
  sheet = "Compound Class"
)
dim(classification_pivot_tbl)
# [1] 106 2

###################
## GENERATE FIGS ##
###################

## fig_1a - Distribution of total abundance & ubiquity

# We will create a data frame which holds the total volatile ubiquity and total
# volatile abundance for plotting.
fig_1a_df <- get_aroma_stats_by_volatles(gcms_pheno_noaid_tbl)
head(fig_1a_df)

# write the fig_1a

# Now, we generate a scatter plot which shows the total volatile ubiquity
# against the total volatile abundance
fig_1a_plot <- fig_1a_df %>%
  mutate(Abundance = as.integer(Abundance / 1e3)) %>%
  ggplot(aes(x = Ubiquity, y = Abundance)) +
  geom_point(aes(alpha = 0.5), size = 5) +
  GLOBAL_THEME +
  theme(
    axis.text.y = element_text(hjust = 0.5)
  ) +
  xlab("Total volatile ubiquity") +
  ylab("Total volatile abundnace (x1000 TIC)") +
  theme(
    legend.position = "none"
  )
ggsave(
  "analyses/data-availability/figs/fig_1a.pdf",
  plot   = fig_1a_plot,
  width  = 10,
  height = 5,
  units  = "in",
  device = cairo_pdf
)


#### fig_1d: Distribution of abundance and the volatiles detected
fig_1de_df <- get_aroma_stats_by_samples(gcms_pheno_noaid_tbl)
fig_1d_plot <- ggplot(
    fig_1de_df,
    aes(
      x = reorder(
        Sample, -Ubiquity
      ),
      y = Ubiquity
    )
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.5), fill = "white", colour = "#1B1716"
  ) +
  geom_hline(yintercept = min(fig_1de_df$Ubiquity), colour = "red") +
  GLOBAL_THEME +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max(fig_1de_df$Ubiquity)),
    breaks = seq(0, max(fig_1de_df$Ubiquity), 10)
  ) +
  xlab("Samples") +
  ylab("Number of volatiles detected")
ggsave(
  "analyses/data-availability/figs/fig_1d.pdf",
  plot   = fig_1d_plot,
  bg     = "white",
  width  = 10,
  height = 5,
  units  = "in",
  device = cairo_pdf
)

fig_1e_plot <- ggplot(
  fig_1de_df, aes(x = reorder(Sample, -Abundance), y = Abundance)
  ) +
  geom_bar(stat = "identity", fill = "white", colour = "#1B1716") +
  geom_hline(yintercept = min(fig_1de_df$Abundance), colour = "red") +
  GLOBAL_THEME +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Samples") +
  ylab("Total volatile abundance (TIC)")
ggsave(
  "analyses/data-availability/figs/fig_1e.pdf",
  plot   = fig_1e_plot,
  bg     = "white",
  width  = 10,
  height = 5,
  units  = "in",
  device = cairo_pdf
)

#### Fig 1B & 1C: Distribution of volatile classes

fig_1b_df <- classification_pivot_tbl %>%
  group_by(Classification) %>%
  summarize(Count = n())

fig_1c_df <- fig_1a_df[, c("Name", "Abundance")] %>%
  inner_join(., classification_pivot_tbl, by = c("Name" = "Compound name")) %>%
  group_by(Classification) %>%
  summarize(TotalAbundance = sum(Abundance))

fig_1b_1c_df <- inner_join(fig_1b_df, fig_1c_df)

fig_1b_plot <- fig_1b_1c_df %>%
  mutate(Classification = fct_reorder(Classification, desc(Count))) %>%
  ggplot(aes(fill = Classification)) +
  geom_bar(
    position = "stack",
    stat = "identity",
    aes(x = 1, y = Count, fill = Classification)
  ) +
  scale_fill_manual(
    guide = guide_legend(reverse = TRUE),
    values = c(
      "#ED6A5A", "#F4F1BB", "#023047", "#9BC1BC", "#5CA4A9", "#E6EBE0",
      "#D4E79E", "#922D50", "#3C1B43", "#501537", "#593959", "#6B4B3E",
      "#55D6BE"
    )
  ) +
  coord_flip() +
  GLOBAL_THEME +
  ylab("Number of volatiles detected") +
  xlab("") +
  theme(
    axis.line.y     = element_blank(),
    axis.text.y     = element_blank(),
    axis.text.x     = element_text(size = 18),
    axis.ticks.y    = element_blank(),
    legend.position = "bottom",
    axis.title.x    = element_text(size = 18),
    legend.title    = element_blank(),
    legend.text     = element_text(size = 11)
  )

ggsave(
  "analyses/data-availability/figs/fig_1b.pdf",
  plot   = fig_1b_plot,
  bg     = "white",
  width  = 20,
  height = 5,
  units  = "in",
  device = cairo_pdf
)

#### fig_1c: Distribution of volatile abundances

fig_1c_plot <- fig_1b_1c_df %>%
  mutate(Classification = fct_reorder(Classification, desc(TotalAbundance))) %>%
  ggplot(aes(fill = Classification, x = 1, y = TotalAbundance)) +
  geom_bar(
    position = "stack", stat = "identity",
    aes(x = 1, fill = Classification)
  ) +
  scale_fill_manual(
    guide = guide_legend(reverse = TRUE),
    values = c(
      "#ED6A5A", "#F4F1BB", "#023047", "#9BC1BC", "#501537", "#E6EBE0",
      "#5CA4A9", "#6B4B3E", "#3C1B43", "#922D50", "#D4E79E", "#593959",
      "#55D6BE"
    )
  ) +
  coord_flip() +
  GLOBAL_THEME +
  ylab("Total volatile abundance (TIC)") +
  xlab("") +
  theme(
    axis.line.y     = element_blank(),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.x     = element_text(size = 18),
    legend.position = "bottom",
    legend.title    = element_blank(),
    axis.title.x    = element_text(size = 18),
    legend.text     = element_text(size = 11)
  )
ggsave(
  "analyses/data-availability/figs/fig_1c.pdf",
  plot   = fig_1c_plot,
  bg     = "white",
  width  = 20,
  height = 5,
  units  = "in",
  device = cairo_pdf
)

#### Arranging fig 1 as a multi-panel fig

fig1_plots <- list()
fig1_plots[[1]] <- fig_1a_plot
fig1_plots[[2]] <- ggarrange(
  fig_1b_plot, fig_1c_plot,
  nrow = 2, ncol = 1,
  common.legend = T,
  labels = c("B", "C")
)
fig1_plots[[3]] <- fig_1d_plot
fig1_plots[[4]] <- fig_1e_plot
fig_1_plot <- ggarrange(
  plotlist = fig1_plots, nrow = 2, ncol = 2, labels = c("A", "", "D", "E")
)
ggsave(
  filename = "analyses/data-availability/figs/fig_1.pdf",
  fig_1_plot,
  dpi      = 600,
  bg       = "white",
  width    = 20,
  height   = 10,
  units    = "in",
  device = cairo_pdf
)

### fig_s1: Distribution of Sample Ubiquity for Volatiles

fig_s1_plot <- ggplot(
  fig_1a_df, aes(x = reorder(Name, -Ubiquity), y = Ubiquity)
) +
  geom_bar(stat = "identity", fill = "white", colour = "black") +
  coord_flip() +
  GLOBAL_THEME +
  theme_avenir(
    panel_x = F, panel_y = F, grid = F, axis = F, axis_col = "black", ticks = T
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(hjust = 0.5, vjust = -0.5),
    axis.title.x = element_text(
      size = 13, hjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0)
    ),
    axis.title.y = element_text(
      size = 13, hjust = 0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)
    ),
  ) +
  ylab("Number of samples") +
  xlab("") +
  scale_y_continuous(expand = c(0, 0))
ggsave(
  "analyses/data-availability/figs/fig_s1.pdf",
  fig_s1_plot,
  width  = 600,
  height = 1000,
  units  = "px",
  dpi    = 75,
  device = cairo_pdf
)
