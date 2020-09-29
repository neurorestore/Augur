#' Superimpose cell type prioritizations onto a dimensionality reduction plot
#'
#' Visualize the global landscape of the perturbation response across a
#' single-cell dataset by superimposing cell type prioritizations (Augur
#' AUCs, or their relative rank within the dataset) for each cell type onto
#' a dimensionality reduction
#'
#' @param augur a set of Augur results obtained from \code{\link{calculate_auc}}
#' @param sc a Seurat, monocle3, or SingleCellExperiment object containing
#'   dimensionality reduction coordinates
#' @param mode whether to plot the raw AUCs or their relative ranks; one of
#'   \code{'default'} (raw AUCs) or \code{'ranks'}
#' @param reduction the type of dimensionality reduction to extract from the
#'   \code{sc} object; defaults to 'umap'
#' @param top_n optionally, the number of top prioritized cell types to label
#'   in the plot
#' @param cell_type_col the column of the \code{meta} data frame, or the
#'   metadata container in the \code{Seurat}/\code{monocle} object, that
#'   contains cell type labels for each cell in the gene-by-cell expression
#'   matrix; defaults to \code{"cell_type"}
#'
#' @return a \code{ggplot} object.
#'
#' @importFrom magrittr %<>%
#' @importFrom dplyr mutate arrange desc slice pull left_join group_by summarise
#'   %>% n
#' @importFrom scales rescale
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na
#' @importFrom stats median
#' @importFrom tidyselect all_of
#' @import ggplot2
#'
#' @export
plot_umap = function(augur, sc, mode = c('default', 'rank'), reduction = 'umap',
                     top_n = 0, cell_type_col = "cell_type") {
  mode = match.arg(mode)
  aucs = augur$AUC

  # extract fill (rank % or AUC) and label
  if (mode == 'rank') {
    aucs %<>%
      mutate(rank = rank(auc),
             rank_pct = rank / n(),
             rank_pct = scales::rescale(rank_pct, c(0, 1)),
             fill = rank_pct)
    legend_name = "Rank (%)"
    label_fun = function(x) x * 100
    breaks = c(0, 1)
    color_labels = c(0, 100)
  } else {
    aucs %<>% mutate(fill = auc)
    legend_name = "AUC"
    breaks = range(aucs$auc)
    color_labels = format(breaks, format = 'f', digits = 2)
  }
  
  # get meta data and UMAP coordinates
  if ("Seurat" %in% class(sc)) {
    # confirm Seurat is installed
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("install \"Seurat\" R package for Augur compatibility with ",
           "input Seurat object", call. = FALSE)
    }
    meta = sc@meta.data %>% as.data.frame()
    red_coord = sc@reductions[[reduction]]@cell.embeddings
  } else if ("monocle3" %in% class(sc)) {
    # confirm monocle3 is installed
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      stop("install \"monocle3\" R package for Augur compatibility with ",
           "input monocle3 object", call. = FALSE)
    }
    meta = monocle3::pData(input) %>%
      droplevels() %>%
      as.data.frame()
    reduction = toupper(reduction)
    red_coord = sc@int_colData@listData$reducedDims[[reduction]]
  } else if ("SingleCellExperiment" %in% class(input)) {
    # confirm SingleCellExperiment is installed
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("install \"SingleCellExperiment\" R package for Augur ",
           "compatibility with input SingleCellExperiment object",
           call. = FALSE)
    }
    meta = SummarizedExperiment::colData(input) %>%
      droplevels() %>%
      as.data.frame()
    reduction = toupper(reduction)
    red_coord = sc@int_colData@listData$reducedDims[[reduction]]
  }
  colnames(red_coord)[1:2] = c('coord_x', 'coord_y')

  # assign fill in Seurat object
  meta$auc = aucs$fill[match(meta[[cell_type_col]], aucs$cell_type)]

  # highlight a handful of top-ranking cell types
  labeled_types = aucs %>%
    arrange(desc(fill)) %>%
    head(top_n) %>%
    pull(cell_type)

  # create plotting structure
  plot_data = red_coord %>%
    as.data.frame() %>%
    rownames_to_column(var = 'barcode') %>%
    left_join(meta %>%
                rownames_to_column(var = 'barcode') %>%
                dplyr::select(barcode, all_of(cell_type_col), auc),
              by = 'barcode')
  colnames(plot_data)[colnames(plot_data) == cell_type_col] = 'cell_type'

  # get labels, and location of the labels
  labels = plot_data %>%
    filter(cell_type %in% labeled_types) %>%
    group_by(cell_type) %>%
    summarise(coord_x = median(coord_x),
              coord_y = median(coord_y),
              auc = mean(auc)) %>%
    drop_na()

  # create plot
  size_sm = 6
  size_lg = 7
  xlab = paste(ifelse(reduction == 'umap', 'UMAP', reduction), 1)
  ylab = paste(ifelse(reduction == 'umap', 'UMAP', reduction), 2)
  p = ggplot(plot_data, aes(x = coord_x, y = coord_y,
                            color = auc, fill = cell_type)) +
    geom_point(size = 0.4, stroke = 0.0, shape = 16) +
    labs(x = xlab, y = ylab) +
    guides(fill = FALSE,
           color = guide_colorbar(nbin = 10, raster = FALSE, ticks = FALSE,
                                  title.position = 'top', title.hjust = 0.5)) +
    geom_text_repel(data = labels,
                    aes(label = cell_type), color = "black", size = 2,
                    segment.size = 0, box.padding = 0.5,
                    min.segment.length = 0.33) +
    scale_color_distiller(palette = "RdGy", name = legend_name,
                          # breaks = c(0, 0.5, 1), limits = c(0, 1),
                          labels = color_labels,
                          breaks = breaks,
                          na.value = 'white') +
    guides(color = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
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
