#' Compare two cell type prioritizations as a scatterplot
#'
#' Compare two sets of cell type prioritization results, calculated for the 
#' same cell types, by comparing them in a scatterplot, with the AUCs
#' from the first set of Augur results on the x-axis and the second set on the
#' y-axis. This function can be useful for understanding how a particular 
#' parameter or preprocessing choice influences cell type prioritization.
#'
#' @param augur1 first set of Augur results, obtained from 
#'   \code{\link{calculate_auc}}
#' @param augur2 second set of Augur results, obtained from 
#'   \code{\link{calculate_auc}}
#' @param top_n optionally, add labels for the top n cell types whose AUCs 
#'   display the greatest change between the first and second Augur results 
#'   into the plot
#'
#' @return a \code{ggplot2} object
#'
#' @importFrom dplyr mutate left_join arrange desc %>%
#' @import ggplot2
#' @import ggrepel
#' @import pals
#'
#' @export
plot_scatterplot = function(augur1, augur2, top_n = 0) {
  auc1 = augur1$AUC
  auc2 = augur2$AUC
  
  # combine
  df = left_join(auc1, auc2, by = 'cell_type') %>%
    mutate(delta = auc.y - auc.x)
  labels = df %>% 
    arrange(desc(abs(delta))) %>%
    head(top_n)
  
  # plot
  size_sm = 6
  size_lg = 7
  pal = pals::coolwarm(n = 100)
  range = range(df$delta)
  limits = c(-max(abs(range)), max(abs(range)))
  limit_labs = format(limits, format = 'f', digits = 2)
  p = df %>%
    ggplot(aes(x = auc.x, y = auc.y, fill = delta)) +
    geom_abline(aes(intercept = 0, slope = 1), linetype = 'dotted',
                size = 0.3) +
    geom_point(size = 0.8, shape = 21, color = 'black', stroke = 0.33) + 
    geom_text_repel(data = labels, aes(label = cell_type), size = 2,
                    segment.size = 0.3, min.segment.length = 0.2, 
                    color = 'black') +
    labs(x = 'AUC 1', y = 'AUC 2') +
    scale_fill_gradientn(expression(Delta~AUC), colours = pal,
                          limits = limits, breaks = limits, 
                          labels = limit_labs) +
    guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_text(size = size_lg),
          panel.grid = element_blank(),
          strip.text = element_text(size = size_lg),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          legend.position = 'right',
          legend.justification = 'bottom',
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.25, "lines"),
          legend.margin = margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5),
          aspect.ratio = 1)
  p
}
