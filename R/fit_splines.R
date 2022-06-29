#' Title
#'
#' @param training_set xxx
#' @param test_set xxx
#' @param nknots xxx
#' @param norder xxx
#'
#' @return xxx
#' @noRd
#' @importFrom fda create.bspline.basis eval.basis
#' @importFrom stats lm coef
fit_splines <- function(training_set, test_set, nknots, norder) {
  breaks <- unique(quantile(training_set[, 2], probs = seq(0, 1, length = nknots)))
  basis <- create.bspline.basis(rangeval = range(breaks),
                                breaks = breaks,
                                norder = norder)
  m <- eval.basis(training_set[, 2], basis)

  D_X <- cbind(m, as.matrix(training_set[,-c(1:2)]))
  D_Y <- training_set[, 1]

  coefs <- coef(lm(D_Y ~ D_X - 1))

  return(pred_splines(object = list(m = m, basis = basis, coefs = coefs), newdata = test_set[, -1]))
}
