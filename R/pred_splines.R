#' Title
#'
#' @param object xxx
#' @param newdata xxx
#'
#' @return xxx
#' @noRd
#' @importFrom fda eval.basis
pred_splines <- function(object, newdata) {
  x <- as.matrix(newdata)
  n <- nrow(x)
  p <- ncol(x)
  lowdiff <- highdiff <- matrix(FALSE, nrow = nrow(newdata), ncol = 1)
  x.cut <- matrix(0, nrow = n, ncol = 1)
  x.current <- x[, 1]
  x.cut[, 1] <- x.current
  bas <- object$basis
  lower.end <- bas$rangeval[1]
  upper.end <- bas$rangeval[2]
  ind.lower <- x.current < lower.end
  ind.upper <- x.current > upper.end

  lowdiff[ind.lower, 1] <- (x.current - lower.end)[ind.lower]
  highdiff[ind.upper, 1] <- (x.current - upper.end)[ind.upper]

  x.cut[ind.lower, 1] <- lower.end
  x.cut[ind.upper, 1] <- upper.end

  ## Get the slopes at the boundaries
  m <- eval.basis(x.cut[, 1], bas)
  deriv.info <- eval.basis(c(lower.end, upper.end), bas, Lfdobj = 1)
  df <- ncol(m)
  lower.slopes <- deriv.info[1, ]
  upper.slopes <- deriv.info[2, ]

  beta <- object$coefs

  pred.pre <- cbind(m, x[, -1]) %*% beta

  ## Put the design matrix of the first derivates (lower.slopes,
  ## upper.slopes) into one long vector each (same length as index) and
  ## multiply with beta vector and take the sum. I.e. perform the matrix
  ## operation in a bit a special form.
  ## The result are the derivatives at the left- and the right-hand side
  ## boundaries of the training range (of the fitted object with the
  ## current coefficients)
  slopes.left <- rowsum(lower.slopes * beta[2:(df + 1)], group = rep(1, df))
  slopes.right <- rowsum(upper.slopes * beta[2:(df + 1)], group = rep(1, df))

  ## Now we have to multiply the derivatives with the difference
  ## in the x-values (contained in lowdiff and highdiff)
  ## lowdiff and highdiff are matrices with dimensions n x p, i.e. the
  ## dimension of the newdata object.
  ## Each column of lowdiff and highdiff is multiplied with the slope
  ## value. The result will be what we have to add beyond the boundaries.
  ## add.left and add.right will also have dimension n x p.

  ## 'as.array' is here to force a warning message if recycling would
  ## take place (see help file of sweep)
  add.left <- sweep(lowdiff, MARGIN = 2, STATS = as.array(slopes.left), FUN = "*")
  add.right <- sweep(highdiff, MARGIN = 2, STATS = as.array(slopes.right), FUN = "*")

  ## Calculate the final prediction:
  ## Take the prediction of the 'cut-down' matrix and add the linear
  ## extrapolation part (add.left + add.right). We have to take the sum
  ## in each row of the linear extrapolation part (add.left + add.right)
  pred.pre <- pred.pre + rowSums(add.left + add.right)
  return(pred.pre)
}
