#' Performs hyperparameter selection by cross-validation.
#'
#' @param df_treatment_A2 a data frame containing the subset of observations used for fitting the treatment model.
#' @param params_grid a data frame containing the combination of hyperparameter values.
#' @param nfolds numeric, the number of folds for the cross validation.
#'
#'
#' @return a list containing the hyperparameter values chosen by cv and the cv mse.
#' @noRd
#'
#' @importFrom xgboost xgb.DMatrix xgb.cv
#'
#' @examples
#' params_grid <- expand.grid(nrounds = 50,
#' eta = 0.3,
#' max_depth = c(2, 6),
#' subsample = c(0.5, 1),
#' colsample_bytree = 1,
#' early_stopping = TRUE)
#' Z <- rnorm(100)
#' X <- rnorm(100)
#' D <- 2 * Z + X + rnorm(100)
#' Data <- data.frame("D" = D, "Z" = Z, "X" = X)
#' cv_fit <- get_l2boost_parameters(df_treatment_A2 = Data, params_grid = params_grid, nfolds = 10)
#' cv_fit$params_A2
get_l2boost_parameters <- function(df_treatment_A2, params_grid, nfolds) {
  n_A2 <- NROW(df_treatment_A2)

  xgbD_A2 <- xgb.DMatrix(as.matrix(df_treatment_A2[, -1]), label = df_treatment_A2[, 1])

  MSE_CV_A2 <- Inf
  params_A2 <- NULL


  for (i in seq_len(NROW(params_grid))) {
    temp_A2 <- xgb.cv(
      params = list(
        eta = params_grid$eta[i],
        max_depth = params_grid$max_depth[i],
        subsample = params_grid$subsample[i],
        colsample_bytree = params_grid$colsample_bytree[i],
        lambda = 0
      ),
      metrics = "rmse",
      nrounds = params_grid$nrounds[i],
      nfold = 5,
      data = xgbD_A2,
      verbose = FALSE
    )
    mse_cv <- temp_A2$evaluation_log$test_rmse_mean^2

    if (params_grid$early_stopping[i] == TRUE) {
      # Search the boosting iteration with the smallest mse and chooses the smallest number of iteration with an mse not much bigger to save computation time when calculating weight matrix if improvement is very slow.
      m <- min(which(mse_cv <= 1.01 * min(mse_cv)))
    } else {
      m <- params_grid$nrounds[i]
    }
    if (mse_cv[m] <= MSE_CV_A2) {
      params_grid$nrounds[i] <- m
      params_A2 <- params_grid[i, ]
      MSE_CV_A2 <- mse_cv[m]
    }
  }
  return(list(
    "params_A2" = params_A2,
    "mse" = MSE_CV_A2
  ))
}
