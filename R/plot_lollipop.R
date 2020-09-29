#' Plot cell type prioritizations as a 'lollipop' plot
#'
#' Plot the complete ranked list of prioritized cell types as a 'lollipop' plot 
#' (similar to a bar chart, except with each bar replaced by a point and line). 
#' In addition, the exact value of the mean AUC, to three decimal places, is
#' printed in the plot for each cell type. For an example, see
#'  **Extended Data Fig. 7b** from 
#' [doi:10.1038/s41587-020-0605-1](https://dx.doi.org/10.1038/s41587-020-0605-1). 
#'
#' @param augur a set of Augur results obtained from \code{\link{calculate_auc}}
#'
#' @return a \code{ggplot2} object
#'
#' @importFrom dplyr mutate %>%
#' @importFrom stats reorder
#' @import ggplot2
#'
#' @export
plot_lollipop = function(augur) {
  aucs = augur$AUC
  size_sm = 6
  size_lg = 7
  range = range(aucs$auc)
  expand = abs(diff(range)) * 0.1
  p = aucs %>%
    # mutate(auc = ifelse(auc < 0.5, 0.5, auc)) %>%
    ggplot(aes(x = reorder(cell_type, auc), y = auc)) +
    geom_hline(aes(yintercept = 0.5), linetype = 'dotted', size = 0.3) +
    geom_point(size = 0.8) +
    geom_text(aes(label = format(auc, digits = 3),
                  y = ifelse(auc < 0.5, 0.5, auc)), size = 2,
              # nudge_x = 0.35, 
              nudge_y = expand,
              hjust = 0.5) +
    geom_segment(aes(xend = cell_type, yend = 0.5)) +
    scale_y_continuous('AUC', limits = c(min(range[1] - expand, 0.5), 
                                         range[2] + expand * 1.5)) +
    coord_flip() +
    theme_bw() + 
    theme(axis.text.x = element_text(size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = size_lg),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5))
  p
}
