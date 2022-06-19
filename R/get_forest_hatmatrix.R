#' Title
#'
#' @param df_treatment_A1 a data frame containing the subset of observations used to get the hat matrix of the treatment model.
#' @param df_treatment_A2 a data frame containing the subset of observations used to train the treatment model.
#' @param params a list containing the hyperparameters of the boosting method.
#'
#' @return a list containing the following the following components: \tabular{ll}{
#' \code{weight} \tab the hat matrix. \cr
#' \tab \cr
#' \code{forest_A2} \tab the fitted random forest model. \cr
#' }
#' @noRd
#' @importFrom ranger ranger
#' @importFrom stats predict
#'
#' @examples
#' Z <- rnorm(100)
#' X <- rnorm(100)
#' D <- 2 * Z + X + rnorm(100)
#' Data <- data.frame("D" = D, "Z" = Z, "X" = X)
#' df_treatment_A1 <- Data[c(1:50), ]
#' df_treatment_A2 <- Data[-c(1:50), ]
#' params_forest <- list("num_trees" = 200, "mtry" = 1, "max_depth" = 0, "min_node_size" = 5)
#' hat_matrix_forest <- get_forest_hatmatrix(df_treatment_A1 = df_treatment_A1, df_treatment_A2 = df_treatment_A2,
#'   params = params_forest)
#' mean((df_treatment_A1$D - hat_matrix_forest$weight %*% df_treatment_A1$D)^2)
get_forest_hatmatrix <- function(df_treatment_A1, df_treatment_A2, params) {
  forest_A2 <- ranger(D ~ .,
                      data = df_treatment_A2,
                      num.trees = params$num_trees,
                      mtry = params$mtry,
                      max.depth = params$max_depth,
                      min.node.size = params$min_node_size,
                      importance = "impurity"
  )
  leaves <- predict(forest_A2, data = df_treatment_A1, type = "terminalNodes")$predictions
  n_A1 <- NROW(leaves)
  num_trees <- NCOL(leaves)
  forest_hatmatrix <- matrix(0, n_A1, n_A1)
  for (j in seq_len(num_trees)) {
    tree_hatmatrix <- get_tree_hatmatrix(leaves = leaves[, j]) # weight matrix for single tree
    if (params$self_predict == FALSE) {
      diag(tree_hatmatrix) <- 0
      tree_hatmatrix <-
        apply(tree_hatmatrix, 1,
              FUN = function(row_weights) {if (sum(row_weights) > 0) {
                return(row_weights / sum(row_weights))} else {
                  return(row_weights)}})
    }
    forest_hatmatrix <- forest_hatmatrix + tree_hatmatrix / num_trees
  }
  #forest_hatmatrix <- Matrix::Matrix(forest_hatmatrix, sparse = T) # sparse matrix to save memory
  return(list(
    "weight" = forest_hatmatrix,
    "model" = forest_A2
  ))
}
