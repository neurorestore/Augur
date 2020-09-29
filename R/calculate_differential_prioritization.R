#' Perform a statistical test for differential prioritization
#'
#' Execute a permutation test to identify cell types with statistically
#' significant differences in AUC between two different rounds of cell type
#' prioritization (for instance, the response to drugs A and B, as compared
#' to a common untreated control).
#'
#' @param augur1 Augur results from condition 1, obtained from
#'   \code{\link{calculate_auc}}
#' @param augur2 Augur results from condition 2, obtained from
#'   \code{\link{calculate_auc}}
#' @param permuted1 permuted Augur results from condition 1, obtained from
#'   \code{\link{calculate_auc}} with the argument \code{augur_mode = "permute"}
#' @param permuted2 permuted Augur results from condition 2, obtained from
#'   \code{\link{calculate_auc}} with the argument \code{augur_mode = "permute"}
#' @param n_subsamples the number of subsamples to pool when calculating the
#'   mean AUC for each permtation; defaults to 50
#' @param n_permutations the total number of mean AUCs to calculate from a
#'   background distribution
#'
#' @return a data frame containing the following columns:
#' \enumerate{
#'   \item \code{cell_type}: the cell types in the input dataste
#'   \item \code{auc.x}: the AUC in condition 1
#'   \item \code{auc.y}: the AUC in condition 2
#'   \item \code{delta_auc}: the difference in AUCs between conditions
#'   \item \code{b} number of times an equally large difference in AUCs was
#'     observed in the permuted data
#'   \item \code{m}: total number of permutations performed
#'   \item \code{z}: the z score of the observed delta-AUC, relative to the
#'     null distribution
#'   \item \code{pval}: the permutation p-value for the observed delta-AUC
#'   \item \code{padj}: the BH-corrected p-value
#' }
#'
#' @importFrom dplyr left_join mutate group_by summarise n select %>%
#' @importFrom tidyr drop_na
#' @importFrom stats p.adjust
#'
#' @export
calculate_differential_prioritization = function(augur1, augur2,
                                                 permuted1, permuted2,
                                                 n_subsamples = 50,
                                                 n_permutations = 1000) {
  obs1 = augur1$AUC
  obs2 = augur2$AUC

  # first, we need to draw mean AUCs for a number of permutations
  permuted_res1 = permuted1$results
  permuted_res2 = permuted2$results
  # calculate interval size
  n_intervals = max(permuted_res1$subsample_idx) / 50
  # first, average across folds
  permuted_auc1 = permuted_res1 %>%
    filter(metric == 'roc_auc') %>%
    group_by(cell_type, subsample_idx) %>%
    summarise(estimate = mean(estimate, na.rm = TRUE)) %>%
    ungroup()
  permuted_auc2 = permuted_res2 %>%
    filter(metric == 'roc_auc') %>%
    group_by(cell_type, subsample_idx) %>%
    summarise(estimate = mean(estimate, na.rm = TRUE)) %>%
    ungroup()
  # now, draw mean AUCs
  draw_mean_aucs = function(permuted_aucs) {
    seq_len(n_permutations) %>%
      map(~ {
        set.seed(.x)
        permuted_aucs %>%
          mutate(bin = cut(subsample_idx, n_intervals, labels = FALSE)) %>%
          # shuffle the bin randomly and take the first one only
          group_by(cell_type) %>%
          mutate(bin = sample(bin)) %>%
          filter(bin == 1) %>%
          mutate(bin = .x) %>%
          group_by(cell_type, bin) %>%
          summarise(
            mean = mean(estimate, na.rm = T),
            sd = sd(estimate, na.rm = T),
            .groups = 'drop_last'
          ) %>%
          ungroup() %>%
          # change bin to permutation_idx
          dplyr::rename(permutation_idx = bin)
      })
  }
  rnd1 = draw_mean_aucs(permuted_auc1) %>% bind_rows()
  rnd2 = draw_mean_aucs(permuted_auc2) %>% bind_rows()

  # calculate observed delta-AUC
  delta = left_join(obs1, obs2, by = 'cell_type') %>%
    mutate(delta_auc = auc.y - auc.x)
  # calculate permuted delta-AUCs
  rnd = left_join(rnd1, rnd2, by = c('cell_type', 'permutation_idx')) %>%
    mutate(delta_rnd = mean.y - mean.x)

  # summarize
  pvals = delta %>%
    left_join(rnd, by = 'cell_type') %>%
    group_by(cell_type) %>%
    summarise(b = sum(delta_rnd >= delta_auc, na.rm = T), m = n(),
              z = unique((delta_auc - mean(delta_rnd, na.rm = T)) /
                           sd(delta_rnd, na.rm = T)),
              pval = min(2 * (b + 1) / (m + 1),
                         2 * (m - b + 1) / (m + 1))) %>%
    ungroup() %>%
    mutate(padj = p.adjust(pval, 'BH'))

  res = delta %>%
    dplyr::select(cell_type, auc.x, auc.y, delta_auc) %>%
    left_join(pvals, by = 'cell_type') %>%
    drop_na(pval)
  return(res)
}
