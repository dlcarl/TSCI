#' @noRd
#' @importFrom ranger ranger
get_forest_parameters <- function(df_treatment_A2, params_grid) {
  # this functions performs hyper-parameter tuning of the random forest parameters.
  forest_A2 <- NULL
  MSE_oob_A2 <- Inf
  params_A2 <- NULL
  ### use oob error to do hyper-parameter tuning
  for (i in seq_len(nrow(params_grid))) {
    temp_A2 <- ranger(D ~ .,
                      data = df_treatment_A2,
                      num.trees = params_grid$num_trees[i],
                      mtry = params_grid$mtry[i][[1]],
                      max.depth = params_grid$max_depth[i],
                      min.node.size = params_grid$min_node_size[i],
                      importance = "impurity"
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
