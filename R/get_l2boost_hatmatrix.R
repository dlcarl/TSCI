#' Calculates the hat matrix of the boosting fit.
#'
#' @param df_treatment_A1 a data frame containing the subset of observations used to get the hat matrix of the treatment model.
#' @param df_treatment_A2 a data frame containing the subset of observations used to train the treatment model.
#' @param params a list containing the hyperparameters of the boosting method.
#'
#' @return a list containing the following the following components: \tabular{ll}{
#' \code{weight} \tab the hat matrix. \cr
#' \tab \cr
#' \code{model} \tab the fitted boosting model. \cr
#' }
#' @noRd
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom stats predict
#' @importFrom Rfast mat.mult
#'
#' @examples
#' Z <- rnorm(100)
#' X <- rnorm(100)
#' D <- 2 * Z + X + rnorm(100)
#' Data <- data.frame("D" = D, "Z" = Z, "X" = X)
#' df_treatment_A1 <- Data[c(1:50), ]
#' df_treatment_A2 <- Data[-c(1:50), ]
#' params_xgboost <- list("nrounds" = 10, "eta" = 0.3, "max_depth" = 2, "subsample" = 1,
#'   "colsample_bytree" = 1)
#' hat_matrix <- get_l2boost_hatmatrix(df_treatment_A1 = df_treatment_A1, df_treatment_A2 = df_treatment_A2,
#'   params = params_xgboost)
#' mean((df_treatment_A1$D - hat_matrix$weight %*% df_treatment_A1$D)^2)
get_l2boost_hatmatrix <- function(df_treatment_A1,
                                   df_treatment_A2,
                                   params) {
  n_A1 <- NROW(df_treatment_A1)
  n_A2 <- NROW(df_treatment_A2)

  xgbD_A1 <- xgb.DMatrix(as.matrix(df_treatment_A1[, -1]))
  xgbD_A2 <- xgb.DMatrix(as.matrix(df_treatment_A2[, -1]), label = df_treatment_A2[, 1])
  l2boost_A2 <- xgb.train(
    nrounds = params$nrounds,
    params = list(
      eta = params$eta,
      max_depth = params$max_depth,
      subsample = params$subsample,
      colsample_bytree = params$colsample_bytree,
      lambda = params$lambda
    ), data = xgbD_A2)
  nodes_A1 <- predict(l2boost_A2, newdata = xgbD_A1, predleaf = TRUE)


  l2boost_hatmatrix <- matrix(0, n_A1, n_A1)
  f_A2 <- rep(0, n_A2)
  for (iter in seq_len(params$nrounds)) {
    leaves_iter_A1 <- nodes_A1[, iter]
    # calculates the hat matrix
    tree_hatmatrix <- get_tree_hatmatrix(leaves_iter_A1)
    l2boost_hatmatrix <- l2boost_hatmatrix + params$eta *
      mat.mult(tree_hatmatrix, diag(n_A1) - l2boost_hatmatrix)
  }
  return(list(
    "weight" = l2boost_hatmatrix,
    "model" = l2boost_A2
  ))
}
