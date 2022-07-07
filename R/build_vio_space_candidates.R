#' Builds violation space candidates.
#'
#' @param Z Z instrument variable with dimension n by 1.
#' @param vio_space a \code{list} with vectors or matrices as elements or a \code{matrix}.
#' If a \code{list} the sequence of violation spaces to test is build by starting with the empty space
#' and obtaining the next larger violation space by iteratively adding the next
#' element of the list to the current violation space.
#' If a \code{matrix} the sequence of violation spaces to test is build by starting with the empty space
#' and obtaining the next larger violation space by iteratively adding the next column of the matrix
#' to the current violation space.
#'
#' @return A list containing:
#' \describe{
#'     \item{\code{vio_space}}{the \code{matrix} of the largest violation space.}
#'     \item{\code{rm_ind}}{a \code{list} containing the indices to remove to obtain the violation spaces to test for (including the null space).}
#'     \item{\code{Q}}{the number of violation spaces (including the null space).}
#' @noRd
#'
#' @importFrom stats poly
#'
build_vio_space_candidates <- function(Z, vio_space) {
  n <- NROW(Z)
  if (class(vio_space)[1] == "list") {
    vio_space <- lapply(vio_space, as.matrix)
    Q <- length(vio_space) + 1
    v_len <- c(1, sapply(vio_space, NCOL))
    # the indices to remove to identify violation space
    rm_ind <-
      lapply(seq_len(Q - 1), FUN = function(i) {
        seq(from = sum(v_len[seq_len(i)]), to = sum(v_len[seq_len(Q - 1) + 1]))
      })
    # merges the list of violation space to a matrix
    vio_space <- Reduce(cbind, vio_space)
  } else if (class(vio_space)[1] == "matrix") {
    Q <- ncol(vio_space) + 1
    rm_ind <- lapply(seq_len(Q - 1), FUN = function(i) seq(from = i, to = Q - 1))
  }
  return(list(
    vio_space = vio_space,
    rm_ind = rm_ind,
    Q = Q
  ))
}