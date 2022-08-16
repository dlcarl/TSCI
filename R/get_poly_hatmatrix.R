#' @noRd
#' @importFrom stats lm coef
get_poly_hatmatrix <- function(df_treatment_A1,
                               params_grid,
                               order_selection_method,
                               max_iter,
                               conv_tol,
                               gcv,
                               nfolds,
                               ...) {
  # this function is used to select the optimal order of the polynomials
  # and to calculate the hat matrix of the treatment model for A1
  # using a polynomial basis expansion.

  # selects optimal order of the polynomials.
  if (max(sapply(seq_len(length(params_grid)), FUN = function(i) length(params_grid[[i]]))) > 1) {
    if (order_selection_method == "backfitting") {
      optimal_order <- backfitting(df = df_treatment_A1,
                                   params_list = params_grid,
                                   max_iter = max_iter,
                                   conv_tol = conv_tol,
                                   gcv = gcv,
                                   nfolds = nfolds)
    }
    else if (order_selection_method == "grid search") {
      optimal_order <- grid_search(df = df_treatment_A1,
                                   params_list = params_grid,
                                   gcv = gcv,
                                   nfolds = nfolds)
    } else {
      stop("No valid order selection method chosen.")
    }
  } else {
    optimal_order <- unlist(params_grid)
  }

  # performs k-fold cross-validation for the chosen order of polynomials to get an out-of-sample MSE.
  A <- matrix(1, nrow = NROW(df_treatment_A1))
  for (j in seq_len(length(optimal_order))) {
    A <- cbind(A, poly(df_treatment_A1[, j + 1], degree = optimal_order[j]))
  }
  MSE_CV <- cross_validate(Y = df_treatment_A1[, 1],
                           X = A,
                           nfolds)

  A <- matrix(1, nrow = NROW(df_treatment_A1))
  for (j in seq_len(length(optimal_order))) {
    A <- cbind(A, poly(df_treatment_A1[, j + 1], degree = optimal_order[j]))
  }
  poly_hatmatrix <- A %*% chol2inv(chol(t(A) %*% A)) %*% t(A)

  return(list(
    weight = poly_hatmatrix,
    mse = MSE_CV
  ))
}
