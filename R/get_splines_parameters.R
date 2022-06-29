#' Title
#'
#' @param df_treatment xxx
#' @param params_grid xxx
#' @param nfolds xxx
#'
#' @return xxx
#' @noRd
get_splines_parameters <- function(df_treatment, params_grid, nfolds) {
  n <- NROW(df_treatment)

  MSE_CV <- Inf
  params <- NULL

  Resample <- sample(n)
  sub_obs <- floor(n / nfolds)
  for (i in seq_len(NROW(params_grid))) {
    mse_tmp <- rep(NA, nfolds)
    for (j in seq_len(nfolds)) {
      if (j == nfolds) {
        sub_test_index <- seq((j - 1) * sub_obs + 1, n)
      } else {
        sub_test_index <- seq((j - 1) * sub_obs + 1, j * sub_obs)
      }

      sub_test_index <- Resample[sub_test_index]

      D_pred <- fit_splines(training_set = df_treatment[-sub_test_index, ],
                               test_set = df_treatment[sub_test_index, ],
                               nknots = params_grid$nknots[i],
                               norder = params_grid$norder[i])

      mse_tmp[j] <- mean((D_pred - df_treatment[sub_test_index, ]$D)^2)

    }


    mse_cv <- mean(mse_tmp)

    if (mse_cv <= MSE_CV) {
      params <- params_grid[i, ]
      MSE_CV <- mse_cv
    }
  }
  return(list(
    "params" = params,
    "mse" = MSE_CV
  ))
}
