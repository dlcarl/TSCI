#' Calculates the p-Values
#'
#' @param Coef an vector of coefficient estimates.
#' @param SE an vector standard errors of the same size as \code{Coef}.
#' @param beta_test the null hypothesis
#'
#' @return the p-Values
#' @noRd
#' @importFrom stats pnorm
p_val <- function(Coef, SE, beta_test) {
  diff_norm <- abs(Coef - beta_test) / SE
  p_val <- 2 * pnorm(diff_norm, lower.tail = FALSE)
  return(p_val)
}
