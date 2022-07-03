#' xxx
#'
#' @param object xxx
#' @param level xxx
#' @param ... xxx
#'
#' @return xxx
#' @export
confint.tsci <- function(object, level = 0.95, ...) {
  stopifnot(inherits(object, "tsci"))
  alpha <- 1 - level
  if (level == (1 - object$alpha)) {
    conf_intervals <- t(cbind(object$CI_robust, object$CI_all))
  } else if (object$mult_split_method == "FWER") {
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
  } else {
    lower <- c(object$Coef_robust, object$Coef_all) - qnorm(1 - alpha / 2) * c(object$sd_robust, object$sd_all)
    upper <- c(object$Coef_robust, object$Coef_all) + qnorm(1 - alpha / 2) * c(object$sd_robust, object$sd_all)
    conf_intervals <- cbind(lower, upper)
  }
  colnames(conf_intervals) <- c(paste(100 * alpha/2, "%"), paste(100*(1 - alpha/2), "%"))
  rownames(conf_intervals) <- names(c(object$Coef_robust, object$Coef_all))

  return(conf_intervals)
}
