#' @noRd
#' @importFrom ranger ranger
#' @importFrom stats predict
get_forest_hatmatrix <- function(df_treatment_A1, df_treatment_A2, params) {
  # this function calculates the hat matrix of the treatment model for A1 using random forest.
  forest_A2 <- ranger(D ~ .,
                      data = df_treatment_A2,
                      num.trees = params$num_trees,
                      mtry = params$mtry[[1]],
                      max.depth = params$max_depth,
                      min.node.size = params$min_node_size,
                      importance = "impurity"
  )
  leaves <- predict(forest_A2, data = df_treatment_A1, type = "terminalNodes")$predictions
  n_A1 <- NROW(leaves)
  forest_hatmatrix <- matrix(0, n_A1, n_A1)
  for (j in seq_len(params$num_trees)) {
    # weight matrix for single tree
    tree_hatmatrix <- get_tree_hatmatrix(leaves = leaves[, j],
                                         self_predict = params$self_predict)
    # updating weight matrix of the tree
    forest_hatmatrix <- forest_hatmatrix + tree_hatmatrix / params$num_trees
  }
  return(list(
    "weight" = forest_hatmatrix,
    "model" = forest_A2
  ))
}
