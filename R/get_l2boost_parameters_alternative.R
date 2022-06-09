#' Performs hyperparameter selection by cross-validation.
#'
#' @param Data_A2 a data frame containing the subset of observations used for fitting the treatment model.
#' @param params_grid a data frame containing the combination of hyperparameter values.
#' @param nfolds numeric, the number of folds for the cross validation.
#'
#'
#' @return a list containing the hyperparameter values chosen by cv and the cv mse.
#' @noRd
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
#' cv_fit <- get_l2boost_parameters(Data_A2 = Data,
#'   params_grid = params_grid,
#'   nfolds = 10)
#' cv_fit$params_A2
get_l2boost_parameters_alternative <- function(Data_A2, params_grid, nfolds) {
  n_A2 <- NROW(Data_A2)
  # Initialize column in data frame containing the residuals at each iteration of the boosting algorithm
  Data_A2$U <- 0
  # Column names of variables we use to predict D
  covariates <- which(!(colnames(Data_A2) %in% c("D", "U")))
  # Position of column containing the residuals at each iteration of the boosting algorithm
  pos_U <- which(colnames(Data_A2) == "U")
  # Randomize sample. Only relevant if we use a supsample of the data for each tree
  observations <- sample(seq_len(n_A2))
  # Number of observation in each fold. Probably needs some more checks to prevent  unreasonable fold sizes (e.g. 0 or n_A2)
  fold_size <- round(n_A2 / nfolds)


  MSE_CV_A2 <- Inf
  params_A2 <- NULL


  for (i in seq_len(NROW(params_grid))) {
    # initialize list for in-sample estimation of f (expected value of D given the covariates and the instrument) at current boosting iteration
    f_cv <- as.list(rep(0, nfolds))
    # initialize lists for out-of-sample predictions, out-of-sample residuals and out-of-sample mean squared error.
    pred_cv <- as.list(rep(0, nfolds))
    res_cv <- vector("list", length = nfolds)
    mse_cv <- rep(NA, params_grid$nrounds[i])
    for (iter in seq_len(params_grid$nrounds[i])) {
      for (k in seq_len(nfolds)) {
        observations_fold <- observations[((k - 1) * fold_size + 1):min(c(k * fold_size, n_A2))]
        # data frame containing observations in the k-th fold. Used to predict out-of-sample
        df_k_fold <- Data_A2[observations_fold, ]
        # data frame containing observations not in the k-th fold. Used to fit model
        df_other_folds <- Data_A2[-1 * observations_fold, ]
        # Residuals at current iteration
        df_other_folds$U <- df_other_folds$D - f_cv[[k]]
        n_cv <- NROW(df_other_folds)
        # creates a subsample for each tree with number of observations and variables specified by the hyperparameters
        df_tmp <- df_other_folds[
          sample(seq_len(n_cv), round(params_grid$subsample[i] * n_cv)),
          c(pos_U, sample(covariates, params_grid$colsample_bytree[i] * length(covariates)))
        ]
        # fits the tree
        temp_A2 <- rpart::rpart(U ~ .,
                                data = df_tmp,
                                control = rpart::rpart.control(maxdepth = params_grid$max_depth[i])
        )
        # in-sample and out-of-sample predictions
        predictions_in_sample <- stats::predict(temp_A2, newdata = df_other_folds, type = "vector")
        predictions_out_sample <- stats::predict(temp_A2, newdata = df_k_fold, type = "vector")
        f_cv[[k]] <- f_cv[[k]] + params_grid$eta[i] * predictions_in_sample
        pred_cv[[k]] <- pred_cv[[k]] + params_grid$eta[i] * predictions_out_sample
        res_cv[[k]] <- df_k_fold$D - pred_cv[[k]]
      }
      mse_cv[iter] <- mean(unlist(res_cv)^2)
    }

    if (params_grid$early_stopping[i] == TRUE) {
      # Search the boosting iteration with the smallest mse and chooses the smallest number of iteration with an mse not much bigger to save computation time when calculating weight matrix if improvement is very slow.
      m <- min(which(mse_cv <= 1.01 * min(mse_cv)))
    } else {
      m <- params_grid$nrounds[i]
    }
    if (mse_cv[m] <= MSE_CV_A2) {
      ind <- i
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
