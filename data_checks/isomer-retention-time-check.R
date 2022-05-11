# Title     : Isomer Retention Time Test  
# Objective : This script checks the clustering of the compounds based on their
#             retention time
# Created by: tayabsoomro
# Created on: 2021-05-14

####################
## 0.1. LIBRARIES ##
####################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

#######################
## 0.2. DATA IMPORTS ##
#######################

# import the cleaned dataset
cleaned.df <- read.table(
  '../data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated.tsv',
  header = TRUE
)

####################
## 0.3. FUNCTIONS ##
####################

clean_names <- function(name_vector) {
  #' clean up the compound names. Particularly, remove the E/Z R/S configurations,
  #' remove unnecessary punctuation and remove all the spaces. Update the name
  #' column in the dataframe
  #' @param name_vector a vector of strings containing all the names
  cleaned_names <- tolower(
    gsub(
      "[[:space:]]","",
      gsub(
        "[[:punct:]]","",
        sub("\\([EZSR,0-9]{1,}\\)","",name_vector)
      )
    )
  )
  return(cleaned_names)
}

gen_plot <- function(dataf, main_title, x_bound, y_bound, ycol, include_count_plot) {
  #' Generates the retention time plot (1st dim. time vs. 2nd dim. time) for the 
  #' data for a given compound.
  #' 
  #' @param dataf the data frame containing First, Second and CName columns
  #' @param main_title the main title of the graph
  #' @param x_bound the minimum and maximum limits of x-axis
  #' @param y_bound the minimum and maximum limits of y-axis
  #' @param ycol the data to use for y-axis. Either Area or Second dimension time
  #' @param include_count_plot Boolean value indicating whether the geom_count 
  #'                           should be added.
  #' 
  
  plt <- NULL
  
  # figure out whether to put Area or Second dimension time on y-axis
  if (ycol == "Area") {
    plt <- ggplot(dataf, aes_string("First","Area", color = "CName"))
      
  } else if (ycol == "CName") {
    plt <- ggplot(dataf, aes_string("First", "Second", color = "CName")) + 
      ylim(y_bound[1], y_bound[2])
  } else {
    stop(paste0("ERROR: ", ycol, " is not a valid column"))
  }
  
  # add the scatter plot
  plt <- plt + geom_point()
  
  # figure out whether to add geom_count
  if (include_count_plot) plt <- plt + geom_count()
  
  plt <- plt + labs(title = main_title) + 
      xlim(x_bound[1], x_bound[2]) + 
      theme_classic() + 
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      ) + 
      annotate("text", 
               x = 750, 
               y = 0.25, 
               label = paste("N=", nrow(dataf)))  
  
  return(plt)
}

##################################
## 1. PREPARE DATA FOR PLOTTING ##
##################################

# add cleaned names to the data frame
cleaned.df$CName <- clean_names(cleaned.df$Name)

# condense the data frame for plotting
for_plot_cleaned.df <- cleaned.df[,c(3,4,15,6)]
colnames(for_plot_cleaned.df) <- c("First","Second","CName", "Area")

# split the data frame by the compound name
X <- split(for_plot_cleaned.df,for_plot_cleaned.df$CName)

# sort the data from the compound with highest number of rows to the lowest.
X_sorted <- X[order(sapply(X, nrow), decreasing = TRUE)]

# gather the first 109 compounds into a final cleaned table
final_109_compound_names <- unique(names(X_sorted))[1:109]
final_109_cmpds <- as.data.frame(matrix(, nrow = 0, ncol = ncol(cleaned.df)))
colnames(final_109_cmpds) <- names(cleaned.df)
for (name in final_109_compound_names) {
  
  final_109_cmpds <- 
    rbind(
      final_109_cmpds, 
      c(cleaned.df[which(cleaned.df$CName == name),])
    )
}

length(unique(final_109_cmpds$CName))

write.table(
  final_109_cmpds,
  '../data/processed/cleaned/cleaned-gcms-data-with-appleid-mapping-sim-gt-850_manual_curated_109.tsv',
  sep = "\t",
  row.names = FALSE
)

##############################
## 2. PLOTTING ##
##############################

# getting the x and y limits for plotting
x_min <- min(for_plot_cleaned.df$First)
x_max <- max(for_plot_cleaned.df$First)

y_min <- min(for_plot_cleaned.df$Second)
y_max <- max(for_plot_cleaned.df$Second)

# generating graphs for individual compounds
main_idx = 1
for (i in seq(1,length(names(X_sorted)), 16)) {
  print(paste0("Generating ", i, " -> ", i + 15))
  # plots for time vs time
  time_plots_simple <- list()
  time_plots_overlap <- list()
  
  # plots for area vs time
  area_plots_simple <- list()
  area_plots_overlap <- list()
  idx <- 1
  for (j in i:(i + 15)) {
    nom <- names(X_sorted)[j]
    time_plots_simple[[idx]] <- gen_plot(X_sorted[[nom]], nom, c(x_min, x_max), c(y_min, y_max),"CName",FALSE)
    time_plots_overlap[[idx]] <- gen_plot(X_sorted[[nom]], nom, c(x_min, x_max), c(y_min, y_max),"CName",TRUE)
    
    area_plots_simple[[idx]] <- gen_plot(X_sorted[[nom]], nom, c(x_min, x_max), c(y_min, y_max),"Area",FALSE)
    area_plots_overlap[[idx]] <- gen_plot(X_sorted[[nom]], nom, c(x_min, x_max), c(y_min, y_max),"Area",TRUE)
    idx <- idx + 1
  }
  time_gs_simple <- ggarrange(plotlist = time_plots_simple)
  time_gs_overlap <- ggarrange(plotlist = time_plots_overlap)
  
  area_gs_simple <- ggarrange(plotlist = area_plots_simple)
  area_gs_overlap <- ggarrange(plotlist = area_plots_overlap)
  
  # save the time vs time plots
  ggsave(
    paste0('figures/retention/time_vs_time/simple/',i,'-',i + 15,'.png'),
    time_gs_simple, width = 20, height = 10, units = "in", device = 'png'
  )
  ggsave(
    paste0('figures/retention/time_vs_time/overplot/',i,'-',i + 15,'.png'),
    time_gs_overlap, width = 20, height = 10, units = "in", device = 'png'
  )
  
  # save the area vs time plots
  ggsave(
    paste0('figures/retention/area_vs_time/simple/',i,'-',i + 15,'.png'),
    area_gs_simple, width = 20, height = 10, units = "in", device = 'png'
  )
  ggsave(
    paste0('figures/retention/area_vs_time/overplot/',i,'-',i + 15,'.png'),
    area_gs_overlap, width = 20, height = 10, units = "in", device = 'png'
  )
  
  main_idx = main_idx + 1
  
}

