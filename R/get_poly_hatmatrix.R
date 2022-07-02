#' Title
#'
#' @param df_treatment_A1 xxx
#' @param df_treatment_A2 xxx
#' @param params xxx
#'
#' @return xxx
#' @noRd
#' @importFrom fda create.bspline.basis eval.basis
get_poly_hatmatrix <- function(df_treatment_A1,
                                  df_treatment_A2,
                                  params) {
  m <- poly(df_treatment_A1[, 2], degree = params$norder, raw = TRUE, simple = TRUE)
  A <- as.matrix(cbind(1, m, df_treatment_A1[, -c(1,2)]))

  poly_hatmatrix <- A %*% chol2inv(chol(t(A) %*% A)) %*% t(A)

  return(list(
    "weight" = poly_hatmatrix,
    "model" = NULL
  ))
}
