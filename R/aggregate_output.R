#' Output Aggregation
#'
#' @param output_list a list with four elements: \code{value}, \code{error}, \code{warning} and \code{message}.
#' @param alpha the significance level.
#' @param Q the number of tested violation spaces.
#' @param mult_split_method method to for inference if multi-splitting is performed. Either 'DML' or 'FWER'.
#'
#' @return a list with the aggregated statistics.
#' @noRd
#' @importFrom stats median
aggregate_output <- function(output_list, alpha, Q,  mult_split_method) {

  # creates list with NAs
  returnList <- tsci_fit_NA_return(Q)


  # extracts coefficient estimates and standard errors of each iteration.
  Coef_all_matrix <-
    matrix(unlist(lapply(output_list, FUN = function(x) x$value$Coef_all), use.names = FALSE),
           ncol = length(returnList$Coef_all), byrow = TRUE)
  sd_all_matrix <-
    matrix(unlist(lapply(output_list, FUN = function(x) x$value$sd_all), use.names = FALSE),
           ncol = length(returnList$sd_all), byrow = TRUE)
  Coef_robust_matrix <-
    matrix(unlist(lapply(output_list, FUN = function(x) x$value$Coef_robust), use.names = FALSE),
                               ncol = length(returnList$Coef_robust), byrow = TRUE)
  sd_robust_matrix <-
    matrix(unlist(lapply(output_list, FUN = function(x) x$value$sd_robust), use.names = FALSE),
           ncol = length(returnList$sd_robust), byrow = TRUE)


  # extracts and aggregates the statistics that do not depend on mult_split_method.
  returnList$Coef_all[] <- apply(Coef_all_matrix, 2, FUN = median)
  returnList$Coef_robust[] <- apply(Coef_robust_matrix, 2, FUN = median)
  returnList$iv_str[] <-
    apply(matrix(unlist(lapply(output_list, FUN = function(x) x$value$iv_str), use.names = FALSE),
                 ncol = length(returnList$iv_str), byrow = TRUE), 2, median)
  returnList$iv_thol[] <-
    apply(matrix(unlist(lapply(output_list, FUN = function(x) x$value$iv_thol), use.names = FALSE),
                 ncol = length(returnList$iv_thol), byrow = TRUE), 2, median)
  returnList$Qmax[] <-
    apply(matrix(unlist(lapply(output_list, FUN = function(x) x$value$Qmax), use.names = FALSE),
                 ncol = length(returnList$Qmax), byrow = TRUE), 2, sum)
  returnList$q_comp[] <-
    apply(matrix(unlist(lapply(output_list, FUN = function(x) x$value$q_comp), use.names = FALSE),
                 ncol = length(returnList$q_comp), byrow = TRUE), 2, sum)
  returnList$q_robust[] <-
    apply(matrix(unlist(lapply(output_list, FUN = function(x) x$value$q_robust), use.names = FALSE),
                 ncol = length(returnList$q_robust), byrow = TRUE), 2, sum)
  returnList$invalidity[] <-
    apply(matrix(unlist(lapply(output_list, FUN = function(x) x$value$invalidity), use.names = FALSE),
                 ncol = length(returnList$invalidity), byrow = TRUE), 2, sum)


  if (mult_split_method == "FWER") {
    stats_all <-
      lapply(seq_len(NCOL(Coef_all_matrix)),
                           FUN = function(j) {get_FWER_CI(Coef = Coef_all_matrix[, j],
                                                          SE = sd_all_matrix[, j],
                                                          level = 1 - alpha)})
    returnList$pval_all[] <- unlist(lapply(stats_all, FUN = function(x) x$p_value))
    returnList$CI_all[] <- matrix(unlist(lapply(stats_all, FUN = function(x) x$CI)), nrow = 2)

    stats_robust <-
      lapply(seq_len(NCOL(Coef_robust_matrix)),
             FUN = function(j) {get_FWER_CI(Coef = Coef_robust_matrix[, j],
                                            SE = sd_robust_matrix[, j],
                                            level = 1 - alpha)})
    returnList$pval_robust[] <- unlist(lapply(stats_robust, FUN = function(x) x$p_value))
    returnList$CI_robust[] <- matrix(unlist(lapply(stats_robust, FUN = function(x) x$CI)), nrow = 2)

    # this method does not provide standard error estimates.
    returnList$sd_all[] <- NA
    returnList$sd_robust[] <- NA

  } else if (mult_split_method == "DML") {
    stats_all <-
      lapply(seq_len(NCOL(Coef_all_matrix)),
             FUN = function(j) {get_DML_CI(Coef = Coef_all_matrix[, j],
                                            SE = sd_all_matrix[, j],
                                            level = 1 - alpha)})
    returnList$sd_all[] <- unlist(lapply(stats_all, FUN = function(x) x$se))
    returnList$pval_all[] <- unlist(lapply(stats_all, FUN = function(x) x$p_value))
    returnList$CI_all[] <- matrix(unlist(lapply(stats_all, FUN = function(x) x$CI)), nrow = 2)

    stats_robust <-
      lapply(seq_len(NCOL(Coef_robust_matrix)),
             FUN = function(j) {get_DML_CI(Coef = Coef_robust_matrix[, j],
                                           SE = sd_robust_matrix[, j],
                                           level = 1 - alpha)})
    returnList$sd_robust[] <- unlist(lapply(stats_robust, FUN = function(x) x$se))
    returnList$pval_robust[] <- unlist(lapply(stats_robust, FUN = function(x) x$p_value))
    returnList$CI_robust[] <- matrix(unlist(lapply(stats_robust, FUN = function(x) x$CI)), nrow = 2)
  } else {
    # if no valid inference method was selected, remove measures of uncertainty.
    returnList$sd_all[] <- rep(NA, times = NCOL(sd_all_matrix))
    returnList$pval_all[] <- rep(NA, times = NCOL(sd_all_matrix))
    returnList$CI_all[] <- matrix(NA, nrow = 2, ncol = NCOL(Coef_all_matrix))
    returnList$sd_robust[] <- rep(NA, times = NCOL(sd_robust_matrix))
    returnList$pval_robust[] <- rep(NA, times = NCOL(sd_robust_matrix))
    returnList$CI_robust[] <- matrix(NA, nrow = 2, ncol = NCOL(Coef_robust_matrix))
  }
  returnList
}
