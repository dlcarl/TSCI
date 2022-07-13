#' @noRd
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom stats predict
#' @importFrom Rfast mat.mult
get_l2boost_hatmatrix <- function(df_treatment_A1,
                                   df_treatment_A2,
                                   params) {
  # this function calculates the hat matrix of the treatment model for A1 using boosting.
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
