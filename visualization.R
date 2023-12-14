library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)

make_upset_plot <- function(
  # Make upset plot
  df = NULL # dataframe of shape (n_samples, n_sets) showing which sets the samples belong to (NA means not belong to)
){
  l <- lapply(df , function(x) rownames(df)[!is.na(x)])
  m <- make_comb_mat(l, mode = "intersect")
  col <- brewer.pal(length(l),"Set3")
  names(col) = names(l)
  UpSet(m, set_order = names(l), comb_order = order(comb_size(m)), left_annotation = upset_left_annotation(m, add_numbers = TRUE, gp = gpar(fill = col, col = "white")))
}