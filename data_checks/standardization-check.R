# Title     : Check Standardization of Phenotype Table
# Objective : This script checks to see if the standardization of the
#             phenotype table makes sense
# Created by: tayabsoomro
# Created on: 2021-04-06

library(ggplot2)
library(ggpubr)

# import the nursery id and standards area values
std_dat <- read.table('../data/processed/standardized/nursery_id_std_area_vals.tsv')

colnames(std_dat)
# [1] "nursery_id"               "X2.hexanone.1.1.1.3.3.d5"
# [3] "Benzaldehyde.d6"          "Ethyl.acetate.d8"

# hexanone vs. benzaldehyde
plot1 <- ggplot(
  std_dat, aes(x = X2.hexanone.1.1.1.3.3.d5, y = Benzaldehyde.d6)) +
  geom_point() +
  xlab("Hexanone") +
  ylab("Benzaldehyde") +
  theme_classic() +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..rr.label..,
                                                  ..p.label.., sep = "~`,`~")))
# hexanone vs. ethyl acetate
plot2 <- ggplot(std_dat, aes(x = X2.hexanone.1.1.1.3.3.d5, y = Ethyl.acetate.d8)) +
  geom_point() +
  xlab("Hexanone") +
  ylab("Ethyl acetate") +
  theme_classic() +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..rr.label..,
                                                  ..p.label.., sep = "~`,`~")))
# benzaldehyde vs. ethyl acetate
plot3 <- ggplot(std_dat, aes(x = Benzaldehyde.d6, y = Ethyl.acetate.d8)) +
  geom_point() +
  xlab("Benzaldehyde") +
  ylab("Ethyl acetate") +
  theme_classic() +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..rr.label..,
                                                  ..p.label.., sep = "~`,`~")))
# generate the combined correlation figure
ggarrange(plot1, plot2, plot3, nrow = 2, ncol = 2)

