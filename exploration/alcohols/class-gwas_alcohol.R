library(readxl)
library(ggplot2)

source('themes/theme_main.R')

class_gcms_tbl <- read_excel('data/processed/sup_tbl_4-final_class_gcms_phenotype_table.xlsx')
dim(class_gcms_tbl)
# [1] 515 14

# create the distribution of cumulative class abundace plot
for (i in 2:ncol(class_gcms_tbl)) {
  name <- colnames(class_gcms_tbl)[i]
  fname <- paste0(
    "figures/class_gwas_exploration/",
    gsub(",", "",
         gsub("/", "-",
              gsub(" ", "_", name)
         )
    ), "_hist.png"
  )
  dat <- as.numeric(unlist(class_gcms_tbl[, i]))
  png(filename = fname)
  hist(
    dat,
    breaks = 100,
    xlab = paste0(name, " abundance"),
    main = paste0("Distribution of ", name, "\n abudance across samples")
  )
  dev.off()
}


# Alcohol seems to be significant at NAC. Explore!
classification_pivot_tbl <- read_excel('data/processed/classification_pivot.xlsx')
all_alcohol_cmpds <- as.character(unlist(classification_pivot_tbl[classification_pivot_tbl$Classification == "Alcohol", "Name"]))

gcms_pheno_tbl <- read_excel('data/processed/sup_tbl_1-final_gcms_phenotype_table_v2.xlsx')
dim(gcms_pheno_tbl)
# [1] 515 107

alcohol_sums_per_cmpd <- colSums(gcms_pheno_tbl[, all_alcohol_cmpds])
alcohol_sums_per_cmpd.df <- data.frame(Name = names(alcohol_sums_per_cmpd), Value = as.numeric(alcohol_sums_per_cmpd))

ggplot(alcohol_sums_per_cmpd.df, aes(x = reorder(Name, -Value), y = Value)) +
  geom_bar(stat = "identity") +
  GLOBAL_THEME +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ylab("Abundance (TIC)") +
  ggtitle("Distribution of abundance values of all Alcohols")
ggsave(
  filename = "figures/class_gwas_exploration/Distribution_of_Alcohol_abundance_per_compound.png",
  plot = last_plot(),
  bg = "white"
)

# the total abundance of every Alchol compound combined, except for 1-Hexanol
sum(alcohol_sums_per_cmpd.df[alcohol_sums_per_cmpd.df$Name != "1-Hexanol", "Value"])
# [1] 8066.095

# the total abundance of 1-Hexanol
sum(alcohol_sums_per_cmpd.df[alcohol_sums_per_cmpd.df$Name == "1-Hexanol", "Value"])
# [1] 8637.869

# The abundance of 1-Hexanol is higher than all the abundances of other Alcohols combined.