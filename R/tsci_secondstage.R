#' @title Violation Space Selection
#' @description Selects violation space. It allows for user specific hat matrices (\code{weight}). Otherwise call \code{tsci_boosting} or \code{tsci_forest}.
#'
#' @param Y outcome with dimension n by 1
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param vio_space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation space selection will be performed according to provided violation space, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case
#' @param A1_ind the indices of samples in A1
#' @param weight n_A1 by n_A1 weight matrix
#' @param intercept logic, to include intercept in the outcome model or not
#' @param str_thol the minimal value of the threshold of IV strength test
#' @param alpha the significance level
#'
#' @return
#' \describe{
#'     \item{\code{Coef_all}}{a series of point estimators of treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{sd_all}}{standard errors of Coef_all}
#'     \item{\code{CI_all}}{confidence intervals for the treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{Coef_robust}}{the point estimators corresponding to the violation space selected by the robust comparison}
#'     \item{\code{sd_robust}}{the standard errors of Coef_robust}
#'     \item{\code{CI_robust}}{confidence intervals for the treatment effect with the violation space selected by the robust comparison}
#'     \item{\code{iv_str}}{IV strength corresponding to different violation spaces}
#'     \item{\code{iv_thol}}{the threshold of IV strength test corresponding to different violation spaces}
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the first violation space and run TSRF only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q_hat}}{the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' }
#' @export
#' @examples
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n)
#' D <- Z + X + rnorm(n)
#' Y <- D + X + rnorm(n) # no violation
#' A1_ind <- seq_len(n)
#' weight <- X %*% chol(chol2inv(t(X) %*% X)) %*% t(X)
#' tsci_secondstage(Y = Y, D = D, Z = Z, X = X, A1_ind = A1_ind, weight = weight)
#'
#' @importFrom stats coef lm qnorm quantile resid rnorm
tsci_secondstage <- function(Y,
                             D,
                             Z,
                             X,
                             vio_space = NULL,
                             A1_ind,
                             weight,
                             intercept = TRUE,
                             str_thol = 10,
                             alpha = 0.05) {
  Y <- as.matrix(Y)
  D <- as.matrix(D)
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  n <- length(Y)
  n_A1 <- length(A1_ind)

  if (nrow(weight) != ncol(weight)) {
    stop("Transformation matrix must be a square matrix")
  }
  if (nrow(weight) != n_A1) {
    stop("The samples to construct transformation matrix must be the same as samples in A1")
  }

  list_vio_space <- check_vio_space(Z, vio_space)
  vio_space <- list_vio_space$vio_space[A1_ind, ]
  rm_ind <- list_vio_space$rm_ind
  Q <- list_vio_space$Q

  Y_A1 <- Y[A1_ind]
  D_A1 <- D[A1_ind]
  X_A1 <- X[A1_ind]

  outputs <- tsci_selection(
    Y_A1 = Y_A1,
    D_A1 = D_A1,
    X_A1 = X_A1,
    vio_space = vio_space,
    rm_ind = rm_ind,
    Q = Q,
    weight = weight,
    intercept = intercept,
    str_thol = str_thol,
    alpha = alpha
  )

  return(outputs)
}
