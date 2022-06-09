#' Calculates the hat matrix of a tree.
#'
#' @param leaves A character or continuous vector specifying the leaf membership.
#' Entries with the same values are assumed to belong in the same leaf.
#'
#' @return A n x n continuous vector. Hat matrix of the tree.
#' @noRd
#'
#' @examples
#' leaves <- c("A", "B", "B", "C", "B", "B", "A")
#' get_tree_hatmatrix(leaves)
get_tree_hatmatrix <- function(leaves) {
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
