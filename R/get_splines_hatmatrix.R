#' Title
#'
#' @param df_treatment_A1 xxx
#' @param df_treatment_A2 xxx
#' @param params xxx
#'
#' @return xxx
#' @noRd
#' @importFrom fda create.bspline.basis eval.basis
get_splines_hatmatrix <- function(df_treatment_A1,
                                  df_treatment_A2,
                                  params) {
  breaks <- unique(quantile(df_treatment_A1[, 2], probs = seq(0, 1, length = params$nknots)))
  basis <- create.bspline.basis(rangeval = range(breaks),
                                breaks = breaks,
                                norder = params$norder)
  m <- eval.basis(df_treatment_A1[, 2], basis)

  A <- cbind(m, as.matrix(df_treatment_A1[,-c(1:2)]))

  bsplines_hatmatrix <- A %*% chol2inv(chol(t(A) %*% A)) %*% t(A)

  return(list(
    "weight" = bsplines_hatmatrix,
    "model" = NULL
  ))
}
