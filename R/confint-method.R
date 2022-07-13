#' Confidence Intervals of Treatment Effect Estimates for TSCI Fits.
#'
#' @param object a object of class 'tsci'.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names.
#' If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... additional argument(s) for methods.
#'
#' @return a matrix containing the confidence intervals.
#' @export
confint.tsci <- function(object, parm = NULL, level = 0.95, ...) {
  stopifnot(inherits(object, "tsci"))
  alpha <- 1 - level
  conf_intervals <- t(cbind(object$CI_robust, object$CI_all))
  if (is.character(parm)) {
    stopifnot(all(parm %in% rownames(conf_intervals)))
    ind <- which(parm %in% rownames(conf_intervals))
  }
  if (is.numeric(parm)) {
    stopifnot(all(parm %in% seq_len(NROW(conf_intervals))))
    ind <- parm
  }
  if (is.null(parm)) {
    ind <- seq_len(NROW(conf_intervals))
  }

  if (level != (1 - object$alpha) & object$mult_split_method == "FWER") {
    if (is.null(object$coef_all_raw) |
        is.null(object$coef_robust_raw) |
        is.null(object$sd_all_raw) |
        is.null(object$sd_robust_raw)) {
      stop("FWER controled confidence intervals cannot be calculated.
             Rerun the tsci fitting function with output_raw = TRUE")
    } else {
      Coef_matrix <- cbind(object$coef_robust_raw, object$coef_all_raw)
      sd_matrix <- cbind(object$sd_robust_raw, object$sd_all_raw)
      conf_intervals <- t(sapply(seq_len(NCOL(Coef_matrix)),
                                 FUN = function(j) {get_FWER_CI(Coef = Coef_matrix[, j],
                                                                SE = sd_matrix[, j],
                                                                level = level)$CI}))
    }
  }
  if (level != (1 - object$alpha) & object$mult_split_method != "FWER") {
    lower <- c(object$Coef_robust, object$Coef_all) - qnorm(1 - alpha / 2) * c(object$sd_robust, object$sd_all)
    upper <- c(object$Coef_robust, object$Coef_all) + qnorm(1 - alpha / 2) * c(object$sd_robust, object$sd_all)
    conf_intervals <- cbind(lower, upper)
  }
  colnames(conf_intervals) <- c(paste(100 * alpha/2, "%"), paste(100*(1 - alpha/2), "%"))
  rownames(conf_intervals) <- names(c(object$Coef_robust, object$Coef_all))

  return(conf_intervals[ind, ])
}
