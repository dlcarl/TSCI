#' @noRd
get_tree_hatmatrix <- function(leaves) {
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
  w <- 1 / (num.samples)
  weight.vec <- ifelse(ind, yes = w, no = 0)
  weight.mat[ind, ] <- matrix(rep(weight.vec, num.samples), num.samples, byrow = T)
  }
  return(weight.mat)
}
