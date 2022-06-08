#' Calculates the hat matrix of the boosting fit.
#'
#' @param Data_A1 a data frame containing the subset of observations used for fitting the outcome model.
#' @param Data_A2 a data frame containing the subset of observations used for fitting the treatment model.
#' @param params a list containing the hyperparameters of the boosting method.
#' @param xgboost logical, specifies whether xgboost should be used to fit the model.
#'
#' @return a list containing the following the following components: \tabular{ll}{
#' \code{weight} \tab the hat matrix. \cr
#' \tab \cr
#' \code{treeboost_A2} \tab the fitted boosting model. \cr
#' }
#' @export
#'
#' @examples
#' Z <- rnorm(100)
#' X <- rnorm(100)
#' D <- 2 * Z + X + rnorm(100)
#' Data <- data.frame("D" = D, "Z" = Z, "X" = X)
#' Data_A1 <- Data[c(1:50), ]
#' Data_A2 <- Data[-c(1:50), ]
#' params_xgboost <- list("nrounds" = 10, "eta" = 0.3, "max_depth" = 2, "subsample" = 1,
#'   "colsample_bytree" = 1)
#' hat_matrix_xgboost <- get_l2boost_hatmatrix(Data_A1 = Data_A1, Data_A2 = Data_A2,
#'   params = params_xgboost, xgboost = TRUE)
#' params_vanilla <- list("nrounds" = 10, "eta" = 0.3, "max_depth" = 2, "subsample" = 1,
#'   "colsample_bytree" = 1)
#' hat_matrix_vanilla <- get_l2boost_hatmatrix(Data_A1 = Data_A1, Data_A2 = Data_A2,
#'   params = params_vanilla, xgboost = FALSE)
#' mean((Data_A1$D - hat_matrix_xgboost$weight %*% Data_A1$D)^2)
#' mean((Data_A1$D - hat_matrix_vanilla$weight %*% Data_A1$D)^2)
get_l2boost_hatmatrix <- function(Data_A1,
                                   Data_A2,
                                   params,
                                   xgboost) {
  n_A1 <- NROW(Data_A1)
  n_A2 <- NROW(Data_A2)

  if (xgboost == TRUE) {
    xgbD_A1 <- xgboost::xgb.DMatrix(as.matrix(Data_A1[, -1]))
    xgbD_A2 <- xgboost::xgb.DMatrix(as.matrix(Data_A2[, -1]), label = Data_A2[, 1])
    l2boost_A2 <- xgboost::xgb.train(
      nrounds = params$nrounds,
      params = list(
        eta = params$eta,
        max_depth = params$max_depth,
        subsample = params$subsample,
        colsample_bytree = params$colsample_bytree,
        lambda = 0
      ), data = xgbD_A2)
    nodes_A1 <- stats::predict(l2boost_A2, newdata = xgbD_A1, predleaf = TRUE)
  } else {
    # Initialize column in data frame containing the residuals at each iteration of the boosting algorithm
    Data_A2$U <- 0
    # Column names of variables we use to predict D
    covariates <- which(!(colnames(Data_A2) %in% c("D", "U")))
    # Position of column containing the residuals we try to predict at each iteration of the boosting alogrithm
    pos_U <- which(colnames(Data_A2) == "U")
    # Randomize sample. Only relevant if we use only a supsample of the data for each tree
    observations <- sample(seq_len(n_A2))
    l2boost_A2 <- params
  }


  l2boost_hatmatrix <- matrix(0, n_A1, n_A1)
  f_A2 <- rep(0, n_A2)
  for (iter in seq_len(params$nrounds)) {
    if (xgboost == TRUE) {
      leaves_iter_A1 <- nodes_A1[, iter]
    } else {
      Data_A2$U <- Data_A2$D - f_A2
      # creates a subsample for each tree with number of observations and variables specified by the hyperparameters
      df_tmp <- Data_A2[
        sample(seq_len(n_A2), round(params$subsample * n_A2)),
        c(pos_U, sample(covariates, params$colsample_bytree * length(covariates)))]
      # fits the tree
      tree_A2 <- rpart::rpart(U ~ ., data = df_tmp, control = rpart::rpart.control(
          maxdepth = params$max_depth))
      predictions_U_A2 <- stats::predict(tree_A2, newdata = Data_A2, type = "vector")
      f_A2 <- f_A2 + params$eta * predictions_U_A2
      leaves_iter_A1 <- stats::predict(tree_A2, newdata = Data_A1, type = "vector")
    }


    # calculates the hat matrix
    tree_hatmatrix <- get_tree_hatmatrix(leaves_iter_A1)
    l2boost_hatmatrix <- l2boost_hatmatrix + params$eta *
      Rfast::mat.mult(tree_hatmatrix, diag(n_A1) - l2boost_hatmatrix)
  }
  return(list(
    "weight" = l2boost_hatmatrix,
    "treeboost_A2" = l2boost_A2
  ))
}
