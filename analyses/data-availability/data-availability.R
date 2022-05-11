# Description: This document contains the code relevant for assessing the various metric
#              such as ubiquity and abundance of the volatile data


#####################
## LIBRARY IMPORTS ##
#####################

library(readxl)
library(tidyverse)
library(grid)
library(ggpubr)
library(forcats)
library("viridis")
source('themes/theme_main.R')


#################
## DEFINITIONS ##
#################

# 1. TotalVolatileUbiquityBySample:
#       This is the number of volatiles that are detected for a given sample.
# 2. TotalVolatileAbundanceBySample:
#       The total cumulative sum of all the volatile concentrations for a given sample
# 3. TotalSampleUbiquityByVolatile:
#       The total number of samples that the given volatile is expressed in

#############################
## DATA LOADING & CURATION ##
#############################

# loading the GCMS phenotype table
gcms_pheno_tbl <- read_excel('data/processed/Supplementary_Data.xlsx', sheet = 'GCMS Data')
dim(gcms_pheno_tbl)
# [1] 515 107

# create a GCMS phenotype table without apple id
gcms_pheno_tbl.noaid <- gcms_pheno_tbl[, 2:ncol(gcms_pheno_tbl)]
dim(gcms_pheno_tbl.noaid)
# [1] 515 106

# load the volatile classification data
classification_pivot <- read_excel('data/processed/Supplementary_Data.xlsx', sheet = 'Compound Class')
dim(classification_pivot)
# [1] 106 2

######################
## GENERATE FIGURES ##
######################

## Figure 1A: Distribution of total abundance & ubiquity

# We will create a data frame which holds the total volatile ubiquity and total volatile abundance for plotting.
TotalSampleUbiquityByVolatiles  <- colSums(gcms_pheno_tbl.noaid != 0)
TotalSampleAbundanceByVolatiles <- colSums(gcms_pheno_tbl.noaid)
fig_1a.df                       <- data.frame(
  Name      = names(TotalSampleUbiquityByVolatiles),
  Ubiquity  = TotalSampleUbiquityByVolatiles,
  Abundance = TotalSampleAbundanceByVolatiles
)
head(fig_1a.df)

# Now, we generate a scatter plot which shows the total volatile ubiquity against the total volatile abundance
fig_1a.plot <- fig_1a.df %>%
  mutate(Abundance = as.integer(Abundance / 1e3)) %>%
  ggplot(aes(x = Ubiquity, y = Abundance)) +
  geom_point(aes(alpha = 0.5), size = 5) +
  GLOBAL_THEME +
  theme(
    axis.text.y = element_text(hjust = 0.5)
  ) + 
  xlab("Total volatile ubiquity") +
  ylab("Total volatile abundnace (TIC)") +
  theme(
    legend.position = "none"
  )
ggsave(
  "figures/data-availability/fig_1a.png",
  plot   = fig_1a.plot,
  bg     = "white",
  width  = 815,
  height = 515,
  dpi    = 75,
  units  = "px"
)

#### Figure 1D: Distribution of volatiles detected
TotalVolatileUbiquityBySample <- rowSums(gcms_pheno_tbl.noaid != 0)
fig_1d.df                     <- data.frame(
  Sample                        = seq_along(TotalVolatileUbiquityBySample),
  TotalVolatileUbiquityBySample = TotalVolatileUbiquityBySample
)
fig_1d.plot                   <-
  ggplot(
    fig_1d.df,
    aes(x = reorder(Sample, -TotalVolatileUbiquityBySample), y = TotalVolatileUbiquityBySample)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5), fill = "white", colour = "#1B1716") +
    geom_hline(yintercept = min(TotalVolatileUbiquityBySample), colour = "red") +
    GLOBAL_THEME +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, max(fig_1d.df$TotalVolatileUbiquityBySample)),
      breaks = seq(0, max(fig_1d.df$TotalVolatileUbiquityBySample), 10)) +
    xlab("Samples") +
    ylab("Number of volatiles detected")
ggsave(
  "figures/data-availability/fig_1d.png",
  plot   = fig_1d.plot,
  bg     = "white",
  width  = 815,
  height = 515,
  dpi    = 75,
  units  = "px"
)

#### Figure 1E: Distribution of volatiles abundance

TotalVolatileAbundnaceBySample <- rowSums(gcms_pheno_tbl.noaid)
fig_1e.df                      <- data.frame(
  Sample    = seq_along(TotalVolatileAbundnaceBySample),
  Abundance = TotalVolatileAbundnaceBySample)
fig_1e.plot                    <- ggplot(fig_1e.df, aes(x = reorder(Sample, -Abundance), y = Abundance)) +
  geom_bar(stat = "identity", fill = "white", colour = "#1B1716") +
  geom_hline(yintercept = min(TotalVolatileAbundnaceBySample), colour = "red") +
  GLOBAL_THEME +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Samples") +
  ylab("Total volatile abundance (TIC)")
ggsave(
  "figures/data-availability/fig_1e.png",
  plot   = fig_1e.plot,
  bg     = "white",
  width  = 815,
  height = 515,
  dpi    = 75,
  units  = "px"
)

#### Figure 1B & 1C: Distribution of volatile classes


fig_1b.df <- classification_pivot %>%
  dplyr::group_by(Classification) %>%
  dplyr::summarize(Count = n())

fig_1c.df <- fig_1a.df[, c("Name", "Abundance")] %>%
  inner_join(., classification_pivot, by = c("Name" = "Compound name")) %>%
  group_by(Classification) %>%
  dplyr::summarize(TotalAbundance = sum(Abundance)) 

fig_1b_1c.df <- inner_join(fig_1b.df, fig_1c.df)

fig_1b.plot <- fig_1b_1c.df %>%
  mutate(Classification = fct_reorder(Classification, desc(Count))) %>%
  ggplot(aes(fill = Classification)) +
  geom_bar(position = "stack", stat = "identity", aes(x = 1, y = Count, fill = Classification)) +
  scale_fill_manual(
    guide  = guide_legend(reverse = TRUE),
    values = c(
      "#ED6A5A", "#F4F1BB", "#023047", "#9BC1BC", "#5CA4A9", "#E6EBE0", "#D4E79E", "#922D50", "#3C1B43", "#501537",
      "#593959", "#6B4B3E", "#55D6BE"
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
  "figures/data-availability/fig_1b.png",
  plot   = fig_1b.plot,
  bg     = "white",
  width  = 2512,
  height = 280,
  dpi    = 75,
  units  = "px"
)


#### Figure 1C: Distribution of volatile abundances

fig_1c.plot <- fig_1b_1c.df %>%
  mutate(Classification = fct_reorder(Classification, desc(TotalAbundance))) %>%
  ggplot(aes(fill = Classification, x = 1, y = TotalAbundance)) +
  geom_bar(position = "stack", stat = "identity", aes(x = 1, fill = Classification)) +
  scale_fill_manual(
    guide  = guide_legend(reverse = TRUE),
    values = c(
      "#ED6A5A", "#F4F1BB", "#023047", "#9BC1BC", "#501537", "#E6EBE0", "#5CA4A9", "#6B4B3E", "#3C1B43", "#922D50", "#D4E79E", "#593959", "#55D6BE"
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
  "figures/data-availability/fig_1c.png",
  plot   = fig_1c.plot,
  bg     = "white",
  width  = 2512,
  height = 280,
  dpi    = 75,
  units  = "px"
)

#### Arranging Figure 1 as a multi-panel figure

fig1_plots      <- list()
fig1_plots[[1]] <- fig_1a.plot
fig1_plots[[2]] <- ggarrange(
  fig_1b.plot, fig_1c.plot,
  nrow          = 2, ncol = 1,
  common.legend = T,
  labels        = c("B", "C")
)
fig1_plots[[3]] <- fig_1d.plot
fig1_plots[[4]] <- fig_1e.plot
figure_1.plot   <- ggarrange(plotlist = fig1_plots, nrow = 2, ncol = 2, labels = c("A", "", "D", "E"))
ggsave(
  filename = "figures/final_figures/Figure_1.svg",
  figure_1.plot,
  dpi      = 600,
  bg       = "white",
  width    = 20,
  height   = 10,
  units    = "in"
)

### Supplementary Figure 1: Distribution of Sample Ubiquity for Volatiles

fig_s1.plot <- ggplot(fig_1a.df, aes(x = reorder(Name, -Ubiquity), y = Ubiquity)) +
  geom_bar(stat = "identity", fill="white", colour="black") +
  coord_flip() +
  GLOBAL_THEME +
  theme_avenir(panel_x = F, panel_y = F, grid = F, axis = F, axis_col = "black", ticks = T) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(hjust = 0.5, vjust = -0.5),
    axis.title.x = element_text(size = 13, hjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 13, hjust = 0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
  ) +
  ylab("Number of samples") +
  xlab("") +
  scale_y_continuous(expand = c(0, 0))
ggsave(
  'figures/final_figures/Figure_S1.svg',
  fig_s1.plot,
  width  = 600,
  height = 1000,
  units  = "px",
  dpi    = 75
)

