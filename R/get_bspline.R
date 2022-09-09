#' @noRd
get_bspline <- function(X, knots = 10) {
  X <- as.matrix(X)
  W <- NULL
  for (i in seq(ncol(X))) {
    X_cur <- X[, i]
    knots_use <- quantile(unique(X_cur), seq(0, 1, length = knots))
    # stopifnot(all.equal(range(knots.use), range(X.cur)))
    basis <- create.bspline.basis(rangeval = range(knots_use), breaks = knots_use, norder = 4)
    basis_val <- eval.basis(X_cur, basis)
    W <- cbind(basis_val, W)
  }
  return(W)
}
