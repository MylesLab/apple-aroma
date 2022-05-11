source("themes/theme_avenir.R")

GLOBAL_PALETTE <- "PuOr"
GLOBAL_DIRECTION <- -1
GLOBAL_ALPHA <- 0.6
GLOBAL_THEME <- theme_avenir(
  panel_x = F,
  panel_y = F,
  grid = F,
  axis = T,
  axis_col = "black",
  ticks = T
) + theme(
  axis.text.x = element_text(size = 10, hjust = 0.5, vjust = -0.5),
  axis.text.y = element_text(size = 10, hjust = -0.5, vjust = 0.5),
  axis.title.x = element_text(size = 13, hjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0)),
  axis.title.y = element_text(size = 13, hjust = 0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
  legend.text = element_text(size = 10)
)
