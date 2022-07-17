#' @noRd
get_tree_hatmatrix <- function(leaves, self_predict = TRUE) {
  # this function returns the hat matrix of a single fitted tree.
  n <- length(leaves)
  # initialize weight matrix for a single tree
  weight.mat <- matrix(0, n, n)
  unique.leaves <- unique(leaves)
  for (i in seq_len(length(unique.leaves))) {
    # indices the observations in the current leaf
    ind <- leaves == unique.leaves[i]
    # the number of observations in the current leaf
    num.samples <- sum(ind)
    if (self_predict) w <- 1 / (num.samples) else w <- 1 / (num.samples - 1)
    weight.vec <- ifelse(ind, yes = w, no = 0)
    weight.mat[ind, ] <- matrix(rep(weight.vec, num.samples), num.samples, byrow = T)
  }
  if (!self_predict) diag(weight.mat) <- 0
  return(weight.mat)
}
