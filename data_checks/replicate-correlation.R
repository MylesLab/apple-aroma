# Title     : Biological Replicates Correlation Test
# Objective : This script looks at the correlation of the two biological
#             replicates for 7 samples.
# Created by: tayabsoomro
# Created on: 2021-03-31

library(ggplot2)
library(ggpubr)

# this is the data frame without connecting the nursery ids to apple ids
cleaned_gcms.df <- read.table('data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated_106.tsv', 1)

# these are the pairs of nursery ids that are biological replicates and we
# would like to see how they correlate with each other.

nids <- list(
  c(6166, 6167),
  c(1008, 1009),
  c(1196, 1197),
  c(2249, 2250),
  c(9304, 9305),
  c(9132, 9134),
  c(8087, 8089)
)

p <- list()
idx <- 1
for (nid in nids) {
  n1 <- cleaned_gcms.df[grep(paste0(".*", nid[1], "_.*"), cleaned_gcms.df$FileName), c('Name', 'Area')]
  n2 <- cleaned_gcms.df[grep(paste0(".*", nid[2], "_.*"), cleaned_gcms.df$FileName), c('Name', 'Area')]

  n1_combined <- aggregate(n1$Area, by = list(Name = n1$Name), FUN = sum)
  n2_combined <- aggregate(n2$Area, by = list(Name = n2$Name), FUN = sum)

  everything <- merge(x = n1_combined, y = n2_combined, by.x = 'Name', by.y = 'Name', all = TRUE)
  colnames(everything) <- c("Name", "Area.x", "Area.y")

  everything$Area.x <- log10(everything$Area.x)
  everything$Area.y <- log10(everything$Area.y)

  diff_df <- everything %>% summarize(Diff = abs(Area.x - Area.y), Name = Name, Area.x = Area.x, Area.y = Area.y)

  p[[idx]] <- ggplot(everything, aes(x = Area.x, y = Area.y)) +
    geom_point() +
    xlab(paste0("Nursery ID: ", nid[1])) +
    ylab(paste0("Nursery ID: ", nid[2])) +
    theme_classic() +
    stat_cor(method = "spearman", aes(label = paste(..r.label.., ..rr.label..,
                                                    ..p.label.., sep = "~`,`~")))
  idx <- idx + 1
  
}
ggarrange(plotlist = p)
