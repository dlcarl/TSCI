#' Confidence Interval and P-Value for Multi-Splitting.
#'
#' @param Coef an vector of coefficient estimates.
#' @param SE an vector standard errors of the same size as \code{Coef}.
#' @param level the confidence level required.
#'
#' @return A list with three elements. \code{se}, \code{CI} and \code{p_value}.
#' @noRd
#' @importFrom stats qnorm pnorm quantile
get_DML_CI <- function(Coef, SE, level) {
  alpha <- 1 - level

  p_val <- function(Coef, SE, beta_test) {
    diff_norm <- abs(Coef - beta_test) / SE
    p_val <- 2 * pnorm(diff_norm, lower.tail = FALSE)
    return(p_val)
  }

  Coef_med <- quantile(Coef, probs = 0.5)
  se <- quantile(sqrt(SE^2 + (Coef - Coef_med)^2), probs = 0.5)
  lower <- Coef_med - qnorm(1 - alpha / 2) * se
  upper <- Coef_med + qnorm(1 - alpha / 2) * se
  p_value <- p_val(Coef_med, se, 0)
  return(list(se = se, CI = c("lower" = lower, "upper" = upper), p_value = p_value))
}
