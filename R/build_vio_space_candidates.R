#' Builds violation space candidates.
#'
#' @param Z Z instrument variable with dimension n by 1.
#' @param vio_space a \code{list} with vectors or matrices as elements.
#' @param create_nested_sequence logical, if \code{TRUE} the sequence of violation spaces to test is build by starting with the empty space
#' and obtaining the next larger violation space by iteratively adding the next
#' element of the \code{vio_space} to the current violation space.
#'
#' @return A list containing:
#' \describe{
#'     \item{\code{vio_space}}{a \code{matrix} consisting of the columns of all violation space candidates.}
#'     \item{\code{vio_ind}}{a \code{list} containing the indices of the violation space candidate.}
#'     \item{\code{Q}}{the number of violation space candidates (including the null space).}
#'     \item{\code{nested_sequence}}{logical. Indicates whether the sequence of violation space candidates is nested.}
#' }
#' @noRd
#'
#' @importFrom stats poly
#'
build_vio_space_candidates <- function(Z, vio_space, create_nested_sequence) {
  n <- NROW(Z)
  vio_space <- lapply(vio_space, as.matrix)
  Q <- length(vio_space) + 1
  # merge violation space candidates into one matrix and remove duplicated columns
  vio_space_matrix <- Reduce(cbind, vio_space)
  duplicated_columns <- duplicated(t(vio_space_matrix))
  vio_space_matrix <- vio_space_matrix[, !duplicated_columns, drop = FALSE]
  # identify the columns of the matrix that belong to each violation space candidate
  vio_ind <- lapply(seq_len(Q - 1), FUN = function(i) {
    which(apply(vio_space_matrix, 2, FUN = function(x) any(colMeans(x == vio_space[[i]]) == 1)))
  })
  vio_space <- vio_space_matrix

  # check if sequence is nested
  nested_sequence <- TRUE
  if (Q > 2) {
    for (i in 2:(Q - 1)) {
      if (!all(vio_ind[[i - 1]] %in% vio_ind[i])) {
        nested_sequence <- FALSE
        break
      }
    }
  }

  if (create_nested_sequence & !nested_sequence) {
    for (i in 2:(Q - 1)) {
      vio_ind[[i]] <- unique(c(vio_ind[[i]], vio_ind[[i - 1]]))
    }
    nested_sequence <- TRUE
  }


  return(list(
    vio_space = vio_space,
    vio_ind = vio_ind,
    Q = Q,
    nested_sequence = nested_sequence
  ))
}

