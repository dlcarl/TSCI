#' hyperparameter tuning using out of bag mean squared error.
#'
#' @param df_treatment_A2 a data frame containing the subset of observations used for fitting the treatment model.
#' @param params_grid a data frame containing the combination of hyperparameter values.
#'
#' @return a list containing the hyperparameter values chosen by cv and the cv mse.
#' @noRd
#'
#' @examples
#' params_grid <- expand.grid(num_trees = 200,
#' mtry = 0.5,
#' max_depth = 0,
#' min_node_size = c(1, 10))
#' Z <- rnorm(100)
#' X <- rnorm(100)
#' D <- 2 * Z + X + rnorm(100)
#' Data <- data.frame("D" = D, "Z" = Z, "X" = X)
#' cv_fit <- get_forest_parameters(df_treatment_A2 = Data, params_grid = params_grid)
#' cv_fit$params_A2
get_forest_parameters <- function(df_treatment_A2, params_grid) {
  forest_A2 <- NULL
  MSE_oob_A2 <- Inf
  params_A2 <- NULL
  ### use oob error to do hyper-parameter tuning
  for (i in seq_len(nrow(params_grid))) {
    temp_A2 <- ranger::ranger(D ~ .,
      data = df_treatment_A2,
      num.trees = params_grid$num_trees[i],
      mtry = params_grid$mtry[i],
      max.depth = params_grid$max_depth[i],
      min.node.size = params_grid$min_node_size[i]
    )
    if (temp_A2$prediction.error <= MSE_oob_A2) {
      forest_A2 <- temp_A2
      params_A2 <- params_grid[i, ]
      MSE_oob_A2 <- temp_A2$prediction.error
    }
  }
  returnList <- list(
    "params_A2" = params_A2,
    "mse" = MSE_oob_A2
  )
  returnList
}
