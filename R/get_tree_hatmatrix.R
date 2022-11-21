#' @noRd
get_tree_hatmatrix <- function(leaves, self_predict = TRUE) {
  # this function returns the hat matrix of a single fitted tree.
  n <- length(leaves)
  # initializes weight matrix for a single tree.
  weight_mat <- matrix(0, n, n)
  unique_leaves <- unique(leaves)
  single_obs <- FALSE
  for (i in seq_len(length(unique_leaves))) {
    # indices the observations in the current leaf.
    ind <- leaves == unique_leaves[i]
    # calculates the number of observations in the current leaf.
    num_samples <- sum(ind)
    # flag if some leaves have only one observation.
    if (num_samples <= 1) single_obs <- TRUE
    # calculates the weights.
    w <- ifelse(self_predict | (num_samples <= 1), 1 / num_samples, 1 / (num_samples - 1))
    weight_vec <- ifelse(ind, yes = w, no = 0)
    weight_mat[ind, ] <- matrix(rep(weight_vec, num_samples), num_samples, byrow = T)
  }
  if (!self_predict) {
    if (single_obs) {
      ind <- rowSums(weight_mat) - diag(weight_mat) > 0
      diag(weight_mat)[ind] <- 0
    } else {
      diag(weight_mat) <- 0
    }
  }
  return(weight_mat)
}
