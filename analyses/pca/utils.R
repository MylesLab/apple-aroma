# Objective : This script contains all the utility functions for PCA

source('themes/theme_avenir.R')

generate_variance_plot <- function(dat, options = NULL) {
    #' Generate the variance plot
    #'
    #' @param dat - data frame with Legend, Variance and Phenotype columns
    #' @param options - plotting options
    #'
    #' @return a ggplot object containing the plot.

  axis_title_size <- options[["axis_title_size"]]
  axis_text_size <- options[["axis_text_size"]]
  legend_text_size <- options[["legend_text_size"]]

  return(
    ggplot(dat, aes(
      fill = Legend, y = Variance, x = reorder(Phenotype, Variance)
    )) +
      geom_bar(
        position = position_stack(reverse = TRUE),
        stat = "identity",
        color = "black"
      ) +
      coord_flip() +
      xlab("") +
      theme_classic2() +
      ylab("% Variance Explained") +
      scale_fill_brewer(palette = "Greens") +
      scale_y_continuous(expand = c(0, 0)) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(
          size=ifelse(is.null(axis_text_size),11,axis_text_size)
        ),
        axis.text.y = element_text(
          size=ifelse(is.null(axis_text_size),11,axis_text_size)
        ),
        axis.title.x = element_text(
          size = ifelse(is.null(axis_title_size),12,axis_title_size)
        ),
        axis.title.y = element_text(
          size = ifelse(is.null(axis_title_size),12,axis_title_size)
        ),
        legend.text = element_text(
          size=ifelse(is.null(legend_text_size), 11, legend_text_size)
        ),
      )
  )
}

generate_pca_biplot <- function(
  dat, choices, color_phenotype, limits, proportion_of_variance, options = NULL
) {
    #' Generate the PCA bi-plot
    #'
    #' @param dat - data frame with Legend, Variance and Phenotype columns
    #' @param choices - a vector representing the PCs for plotting
    #' @param color_phenotype - a vector containing the phenotype that the data
    #'  should be colored by as well as the legend.
    #'  For example c("HarvestDate","Harvest Date (Julian days)")
    #' @param proportion_of_variance - dataframe containing the proportion of
    #'  variances
    #' @param options - any extra options
    #'
    #' @return a ggplot object containing the plot.
  X <- choices[1]
  Y <- choices[2]

  pheno_name <- color_phenotype[1]
  pheno_title <- color_phenotype[2]

  # other parameters
  dot_size <- options[["dot_size"]]
  axis_text_size <- options[["axis_text_size"]]
  axis_title_size <- options[["axis_title_size"]]
  legend_text_size <- options[["legend_text_size"]]
  legend_title_size <- options[["legend_title_size"]]



  ggplot(dat, aes_string(x = X, y = Y)) +
    geom_point(
      aes_string(
        colour = pheno_name
      ),
      size = ifelse(is.null(dot_size),2,dot_size)
    ) +
    scale_color_viridis(
      name = pheno_title,
      option = "magma",
      direction = -1
    ) +
    theme_classic2() + 
    theme(
      legend.position = c(0.9,0.3),
      axis.text.x = element_text(
        size=ifelse(is.null(axis_text_size),11,axis_text_size)
      ),
      axis.text.y = element_text(
        size=ifelse(is.null(axis_text_size),11,axis_text_size)
      ),
      axis.title.x = element_text(
        size=ifelse(is.null(axis_title_size),11,axis_title_size)
      ),
      axis.title.y = element_text(
        size=ifelse(is.null(axis_title_size),11,axis_title_size)
      ),
      legend.title = element_text(
        size=ifelse(is.null(legend_title_size),11,legend_title_size)
      ),
      legend.text = element_text(
        size=ifelse(is.null(legend_text_size), 11, legend_text_size)
      )

    ) +
    guides(colour = guide_colourbar(barwidth = 1, barheight = 10)) +
    xlim(limits) +
    ylim(limits) +
    xlab(sprintf("PC1 (%s%%)", round(as.numeric(proportion_of_variance[1]), 1))) +
    ylab(sprintf("PC2 (%s%%)", round(as.numeric(proportion_of_variance[2]), 1)))
}
