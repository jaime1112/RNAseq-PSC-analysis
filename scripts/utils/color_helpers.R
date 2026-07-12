# Shared color helpers
# Small utilities used by more than one plotting function.

#' Build a named categorical color vector for heatmap annotations
#'
#' Assigns one distinct color per unique value, drawing from a pool of
#' ColorBrewer "Set2" (8) + "Set3" (12) = up to 20 qualitative colors. Returns
#' NULL when there are more levels than the pool can cover, so the caller can
#' fall back to the plotting package's automatic colors.
#'
#' brewer.pal() has a minimum request of 3 and fixed maxima per palette, so the
#' pool is built once at full size and then sliced. This avoids the off-by-count
#' warnings/mis-assignments that arise from requesting fewer than 3 colors, or
#' 8 + a small remainder.
#'
#' @param values Vector of category values; levels are taken as unique(values),
#'   preserving first-appearance order.
#' @return Named character vector of hex colors (names = the unique values), or
#'   NULL if there are more unique values than available colors.
categorical_annotation_colors <- function(values) {
  levels_present <- unique(values)
  n <- length(levels_present)

  color_pool <- c(RColorBrewer::brewer.pal(8, "Set2"),
                  RColorBrewer::brewer.pal(12, "Set3"))

  if (n < 1 || n > length(color_pool)) {
    return(NULL)
  }

  stats::setNames(color_pool[seq_len(n)], levels_present)
}
