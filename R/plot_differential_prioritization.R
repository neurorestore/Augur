#' Plot the results of a differential prioritization analysis
#' 
#' After performing a statistical test for differential prioritization using
#' \code{\link{calculate_differential_prioritization}}, plot the results
#' in a scatterplot, highlighting cell types with significant differences
#' between conditions.
#'
#' @param results the output from 
#'   \code{\link{calculate_differential_prioritization}}
#' @param top_n optionally, the number of top prioritized cell types to label
#'   in the plot
#'
#' @return a \code{ggplot2} object
#' 
#' @importFrom dplyr filter arrange mutate %>%
#' @importFrom tidyr drop_na
#' @import ggplot2
#' @import pals
#' @import ggrepel
#'
#' @export
plot_differential_prioritization = function(results, top_n = 0) {
  labels = results %>%
    filter(pval < 0.05) %>%
    arrange(pval, desc(abs(z))) %>%
    head(top_n)
  pal = pals::kovesi.diverging_cwm_80_100_c22(n = 20)[c(1, 20)] %>% c("grey80")
  size_sm = 6
  size_lg = 7
  p = results %>%
    mutate(color = factor(ifelse(pval < 0.05,
                                 ifelse(z > 0, 'condition 2', 'condition 1'),
                                 'n.s.'),
                          levels = c('condition 1', 'condition 2', 'n.s.'))) %>%
    drop_na(auc.x, auc.y, color) %>%
    ggplot(aes(x = auc.x, y = auc.y, color = color)) + 
    geom_abline(aes(intercept = 0, slope = 1), linetype = 'dotted') +
    geom_point(size = 0.7) + 
    geom_text_repel(data = labels,
                    aes(label = cell_type), color = "black", size = 2,
                    segment.size = 0.2, box.padding = 0.5,
                    min.segment.length = 0.33) +
    scale_color_manual('', values = pal, drop = F) + 
    scale_x_continuous(name = 'AUC 1') +
    scale_y_continuous(name = 'AUC 2') +
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
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5),
          aspect.ratio = 1)
  p
}