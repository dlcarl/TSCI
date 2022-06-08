#' Title
#'
#' @param Data_A1 a data frame containing the subset of observations used for fitting the outcome model.
#' @param Data_A2 a data frame containing the subset of observations used for fitting the treatment model.
#' @param params a list containing the hyperparameters of the boosting method.
#'
#' @return a list containing the following the following components: \tabular{ll}{
#' \code{weight} \tab the hat matrix. \cr
#' \tab \cr
#' \code{forest_A2} \tab the fitted random forest model. \cr
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
#' params_forest <- list("num_trees" = 200, "mtry" = 1, "max_depth" = 0, "min_node_size" = 5)
#' hat_matrix_forest <- get_forest_hatmatrix(Data_A1 = Data_A1, Data_A2 = Data_A2,
#'   params = params_forest)
#' mean((Data_A1$D - hat_matrix_forest$weight %*% Data_A1$D)^2)
get_forest_hatmatrix <- function(Data_A1, Data_A2, params) {
  forest_A2 <- ranger::ranger(D ~ .,
    data = Data_A2,
    num.trees = params$num_trees,
    mtry = params$mtry,
    max.depth = params$max_depth,
    min.node.size = params$min_node_size
  )
  leaves <- stats::predict(forest_A2, data = Data_A1, type = "terminalNodes")$predictions
  n_A1 <- NROW(leaves)
  num_trees <- NCOL(leaves)
  forest_hatmatrix <- matrix(0, n_A1, n_A1)
  for (j in seq_len(num_trees)) {
    tree_hatmatrix <- matrix(0, n_A1, n_A1) # weight matrix for single tree
    unique_leaves <- unique(leaves[, j])
    for (i in seq_len(length(unique_leaves))) {
      ind <- leaves[, j] == unique_leaves[i] # indices of samples in the node
      num_samples <- sum(ind) # number of samples in the node
      w <- 1 / num_samples
      weight_vec <- ifelse(ind, yes = w, no = 0)
      tree_hatmatrix[ind, ] <- matrix(rep(weight_vec, num_samples), num_samples, byrow = T) / num_trees
    }
    forest_hatmatrix <- forest_hatmatrix + tree_hatmatrix
  }
  #forest_hatmatrix <- Matrix::Matrix(forest_hatmatrix, sparse = T) # sparse matrix to save memory
  return(list(
    "weight" = forest_hatmatrix,
    "forest_A2" = forest_A2
  ))
}
