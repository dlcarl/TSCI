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
  params <- unlist(params)
  A <- matrix(1, nrow = NROW(df_treatment_A1))
  for (j in seq_len(length(params))) {
    A <- cbind(A, poly(df_treatment_A1[, j + 1], degree = params[j]))
  }
  poly_hatmatrix <- A %*% chol2inv(chol(t(A) %*% A)) %*% t(A)

  return(list(
    "weight" = poly_hatmatrix,
    "model" = NULL
  ))
}
