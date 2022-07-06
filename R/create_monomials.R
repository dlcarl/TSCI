#' Monomials as Violation Space Candidates
#'
#' @param Z observations of the instrumental variable(s). Either a numeric vector of length n
#' or a numeric matrix with dimension n by s.
#' @param degree The degree up to which monomials should be created.
#' @param type One out of \code{monomials_main} or \code{monomials_full}. \cr
#' \code{monomials_main} creates the monomials for the polynomials of each instrumental variable up to degree \code{degree}. \cr
#' \code{monomials_full} creates the monomials for the polynomials of a combination of all instrumental variables up to degree \code{degree}.
#'
#' @return A ordered list. Each element is a matrix consisting of the monomials to be added to
#' the next violation space candidate.
#'
#' @details assuming there are 3 instrumental variables Z1, Z2, and Z3 and \code{degree} = d,
#' \code{monomials_main} creates the monomials of the polynomials (Z1 + 1)^d, (Z2 + 1)^d, (Z3 + 1)^d without the constants and
#' \code{monomials_full} creates the monomials of the polynomials (Z1 + Z2 + Z3)^1, (Z1 + Z2 + Z3)^2, ..., (Z1 + Z2 + Z3)^d without the constants.
#' Thus \code{type} = \code{monomials_main} does not include interactions between the instrumental variables.
#'
#' @export
#'
#' @examples
#' Z <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
#' vio_space <- create_monomials(Z = Z, degree = 4, type = "monomials_full")
create_monomials <- function(Z, degree, type) {

  # check if input is in the correct format and valid
  if (!is.numeric(Z))
    error_message <- paste(error_message, "Z is not numeric.", sep = "\n")
  if (!is.numeric(degree))
    error_message <- paste(error_message, "degree is not numeric.", sep = "\n")
  if (!(type %in% c("monomials_main", "monomials_full")))
    error_message <- paste(error_message, "No valid type
                           selected. Choose either 'monomials_main' or 'monomials_full'.", sep = "\n")

  # build monomials
  Z <- as.matrix(Z)
  n <- NROW(Z)
  p <- NCOL(Z)
  monomials <- vector("list", length = degree + 1)
  monomials[[1]] <- matrix(1, nrow = n)
  for (q in seq_len(degree)) {
    if (type == "monomials_full") {
      current_monomial <- matrix(NA, nrow = n, ncol = p * NCOL(monomials[[q]]))
      for (i in seq_len(NCOL(monomials[[q]]))) {
        for (j in seq_len(p)) {
          current_monomial[, (i - 1) * p + j] <- monomials[[q]][, i] * Z[, j]
        }
      }
    } else if (type == "monomials_main") {
      if (q == 1) current_monomial <- Z
      else  current_monomial <- monomials[[q]] * Z
    }

    duplicated_columns <- duplicated(t(round(current_monomial, 10)))
    monomials[[q + 1]] <- current_monomial[, !duplicated_columns]
  }
  return(monomials[-1])
}
