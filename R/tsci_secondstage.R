#' Two Stage Curvature Identification with User Provided Hat Matrix
#' @description \code{tsci_secondstage} implements Two Stage Curvature Identification
#' (Guo and Bühlmann 2022) for a user-provided hat matrix. Through a data-dependent way it
#' tests for the smallest sufficiently large violation space among a pre-specified
#' sequence of nested violation space candidates. Point and uncertainty estimates
#' of the treatment effect for all violation space candidates including the
#' selected violation space will be returned amongst other relevant statistics.
#'
#' @param Y observations of the outcome variable. Either a numeric vector of length n
#' or a numeric matrix with dimension n by 1.
#' @param D observations of the treatment variable. Either a numeric vector of length n
#' or a numeric matrix with dimension n by 1.
#' @param Z observations of the instrumental variable(s). Either a numeric vector of length n
#' or a numeric matrix with dimension n by s.
#' @param X observations of baseline covariate(s). Either a numeric vector of length n
#' or a numeric matrix with dimension n by p or \code{NULL}
#' (if no covariates should be included).
#' @param vio_space either a numeric matrix with dimension n by q or a list with
#' numeric vectors of length n and/or numeric matrices with n rows as elements to
#' specify the violation space candidates or \code{NULL}.
#' If a matrix or a list, then the violation space candidates (in form of matrices)
#' are defined sequentially starting with an empty violation matrix and subsequently
#' adding the next column of the matrix or element of the list to the current violation matrix.
#' See Details for more information.
#' @param weight the hat matrix of the treatment model.
#' @param A1_ind indices of the observations that wil be used to fit the outcome model.
#' Must be of same length as the number of rows and columns of \code{weight}.
#' If \code{NULL} all observations will be used.
#' @param intercept logical. If \code{TRUE} an intercept is included in the outcome model.
#' @param str_thol minimal value of the threshold of IV strength test.
#' @param alpha the significance level.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'     \item{\code{Coef_all}}{a series of point estimates of the treatment effect
#'     for the different violation space candidates and the OLS estimate.}
#'     \item{\code{sd_all}}{standard errors of Coef_all.}
#'     \item{\code{pvall_all}}{p-values of the treatment effect estimates for the
#'     different violation space candidates and for the OLS estimate.}
#'     \item{\code{CI_all}}{confidence intervals for the treatment effect for the
#'     different violation space candidates and for the OLS estimate.}
#'     \item{\code{Coef_robust}}{the point estimators of the treatment effect for
#'     the selected violation spaces.}
#'     \item{\code{sd_robust}}{the standard errors of Coef_robust.}
#'     \item{\code{pvall_all}}{p-values of the treatment effect estimates for the
#'     selected violation spaces.}
#'     \item{\code{CI_robust}}{confidence intervals for the treatment effect for
#'     the selected violation spaces.}
#'     \item{\code{iv_str}}{IV strength for the different violation space candidates.}
#'     \item{\code{iv_thol}}{the threshold for the different violation space candidates.}
#'     \item{\code{Qmax}}{a named vector containing the number of times the
#'     violation space candidates were the largest acceptable violation space by
#'     the IV strength test.
#'     The value of the element named "OLS" is the number of times the instrument
#'     was too weak even for the empty violation space.}
#'     \item{\code{q_comp}}{a named vector containing the number of times the
#'     violation space candidates were the selected violation space.
#'     The value of the element named "OLS" is the number of times the instrument
#'     was too weak even for the empty violation space.}
#'     \item{\code{q_robust}}{a named vector containing the number of times the
#'     violation space candidates were the selected violation space by the robust method.
#'     The value of the element named "OLS" is the number of times the instrument
#'     was too weak even for the empty violation space.}
#'     \item{\code{invalidity}}{a named vector containing the number of times
#'     the instrument was considered valid and invalid. The instrument is considered
#'     invalid if the selected violation space is larger than the empty space.}
#'     \item{\code{mse}}{the out-of-sample mean squared error of the treatment model.}
#'}
#'
#' @details The treatment and outcome models are assumed to be of the following forms:
#' \deqn{D_i = g(Z_i, X_i) + \delta_i}
#' \deqn{Y_i = \beta * D_i + h(Z_i, X_i) + \epsilon_i}
#' where \eqn{g(Z_i, X_i)} is assumed to be estimated by the user and
#' \eqn{h(Z_i X_i)} is approximated using the violation space candidates and by
#' a linear combination of baseline covariates. The errors are allowed to be heteroscedastic.
#' \eqn{A1} is used to fit the outcome model. \cr \cr
#' The violation space candidates are required to be in a nested sequence. The specification
#' of suitable violation space candidates is a crucial step because a poor approximation
#' of \eqn{h(Z_i, X_i)} might not address the bias caused by the violation of the IV assumption sufficiently.
#' The function \code{create_monomials} can be used to create such a nested sequence for a
#' predefined type of violation space candidates (monomials). \cr \cr
#'
#' @references
#' \itemize{
#' \item{Zijian Guo, and Peter Bühlmann. Two Stage Curvature Identification with
#' Machine Learning: Causal Inference with Possibly Invalid Instrumental Variables.
#' \emph{arXiv:2203.12808}, 2022}
#' \item{Nicolai Meinshausen, Lukas Meier, and Peter Bühlmann. P-values for high-dimensional
#' regression. \emph{Journal of the American Statistical Association},
#' 104(488):1671-1681, 2009. 16, 18}
#' \item{Victor Chernozhukov, Denis Chetverikov, Mert Demirer, Esther Duflo, Christian Hansen,
#' Whitney Newey, and James Robins. Double/debiased machine learning for treatment
#' and structural parameters: Double/debiased machine learning.
#' \emph{The Econometrics Journal}, 21(1), 2018. 4, 16, 18}
#' }
#' @export
#'
#' @examples
#' # dimension
#' p <- 10
#' # sample size
#' n <- 1000
#' # interaction value
#' inter_val <- 1
#' # the IV strength
#' a <- 1
#' # violation strength
#' tau <- 1
#' f <- function(x) {a * (1 * sin(2 * pi * x) + 1.5 * cos(2 * pi * x))}
#' rho <- 0.5
#' Cov <- stats::toeplitz(rho^c(0 : p))
#' mu <- rep(0, p + 1)
#' # true effect
#' beta <- 1
#' alpha <- as.matrix(rep(-0.3, p))
#' gamma <- as.matrix(rep(0.2, p))
#' inter <- as.matrix(c(rep(inter_val, 5),rep(0, p - 5)))
#'
#'
#' # generate the data
#' mu_error <- rep(0,2)
#' Cov_error <- matrix(c(1, 0.5, 0.5,1), 2, 2)
#' Error <- MASS::mvrnorm(n, mu_error, Cov_error)
#' W_original <- MASS::mvrnorm(n, mu, Cov)
#' W <- pnorm(W_original)
#' # instrument variable
#' Z <- W[, 1]
#' # baseline covariates
#' X <- W[, -1]
#' # generate the treatment variable D
#' D <- f(Z) + X %*% alpha + Z * X %*% inter + Error[, 1]
#' # generate the outcome variable Y
#' Y <- D * beta + tau * Z + X %*% gamma + Error[, 2]
#'
#' # get hat matrix of outcome model
#' A <- cbind(Z, Z^2, Z^3, Z^4, Z*X, X)
#' weight <- A %*% chol2inv(chol(t(A) %*% A)) %*% t(A)
#'
#' # Two Stage L2 Boosting
#' vio_space <- create_monomials(Z, 4, "monomials_main")
#' output <- tsci_secondstage(Y, D, Z, X, vio_space, weight)
#' # point estimates
#' output$Coef_robust
#' # standard errors
#' output$sd_robust
#' # confidence intervals
#' output$CI_robust
#'
#' @importFrom stats coef lm qnorm quantile resid rnorm
tsci_secondstage <- function(Y,
                             D,
                             Z,
                             X,
                             vio_space,
                             weight,
                             A1_ind = NULL,
                             intercept = TRUE,
                             str_thol = 10,
                             alpha = 0.05) {
  # check that input is in the correct format
  error_message <- NULL
  if (!is.numeric(Y))
    error_message <- paste(error_message, "Y is not numeric.", sep = "\n")
  if (!is.numeric(D))
    error_message <- paste(error_message, "D is not numeric.", sep = "\n")
  if (!is.numeric(Z))
    error_message <- paste(error_message, "Z is not numeric.", sep = "\n")
  if (!is.numeric(X) & !is.null(X))
    error_message <- paste(error_message, "X is beither numeric nor NULL.", sep = "\n")
  if (!is.numeric(A1_ind) & !is.null(A1_ind))
    error_message <- paste(error_message, "A1_ind is neither numeric nor NULL.", sep = "\n")
  if (!is.numeric(weight))
    error_message <- paste(error_message, "weight is not numeric.", sep = "\n")
  if (!is.logical(intercept))
    error_message <- paste(error_message, "intercept is neither TRUE nor FALSE.", sep = "\n")
  if (!is.matrix(vio_space) & !is.list(vio_space) & !is.null(vio_space))
    error_message <- paste(error_message, "vio_space is neither a matrix nor a list nor NULL", sep = "\n")
  if (is.matrix(vio_space)) {
    if (!is.numeric(vio_space))
      error_message <- paste(error_message, "vio_space is not numeric", sep = "\n")
  } else if (is.list(vio_space)) {
    if (!is.numeric(unlist(vio_space)))
      error_message <- paste(error_message, "vio_space is not numeric", sep = "\n")
  }
  if (!is.numeric(str_thol))
    error_message <- paste(error_message, "str_thol is not numeric.", sep = "\n")
  if (!is.numeric(alpha))
    error_message <- paste(error_message, "alpha is not numeric.", sep = "\n")

  if (!is.null(error_message))
    stop(error_message)

  # check if inputs are possible
  p <- NCOL(Z) + ifelse(is.null(X), 0, NCOL(X))
  if (length(unique(sapply(list(Y, D, Z), FUN = function(variable) NROW(variable)))) > 1)
    error_message <- paste(error_message, "Y, D and Z have not the same amount of observations.", sep = "\n")
  else {
    n <- NROW(Y)
    if (!is.null(X))
      if(NROW(X) != n)
        error_message <- paste(error_message, "X has not the same amount of observations as Y.", sep = "\n")
    if (is.matrix(vio_space))
      if (NROW(vio_space) != n)
        error_message <- paste(error_message, "vio_space has not the same amount of observations as Y.", sep = "\n")
    if (is.list(vio_space))
      if (length(unique(sapply(vio_space, FUN = function(variable) NROW(variable)))) > 1)
        error_message <- paste(error_message, "vio_space has not the same amount of observations as Y.", sep = "\n")
  }
  if (any(is.na(Y)))
    error_message <- paste(error_message, "There are NA's in Y.", sep = "\n")
  if (any(is.na(D)))
    error_message <- paste(error_message, "There are NA's in D.", sep = "\n")
  if (any(is.na(Z)))
    error_message <- paste(error_message, "There are NA's in Z.", sep = "\n")
  if (!is.null(X))
    if(any(is.na(X)))
      error_message <- paste(error_message, "There are NA's in X.", sep = "\n")
  if (!is.null(A1_ind)) {
    if (any(duplicated(A1_ind)))
      error_message <- paste(error_message, "There are duplicates in A1_ind.", sep = "\n")
    if (length(A1_ind) > NROW(Y))
      error_message <- paste(error_message, "A1_ind contains more elements than Y has observations", sep = "\n")
    if (any(A1_ind < 0) | any(A1_ind > NROW(Y)))
      error_message <- paste(error_message, "A1_ind contains elements that are negative or larger than Y has observations", sep = "\n")
    if (length(A1_ind) != NROW(weight))
      error_message <- paste(error_message, "A1_ind contains not the same amount of elements as weight has rows.", sep = "\n")
  }
  if (NROW(weight) != NROW(Y))
    error_message <- paste(error_message, "weight has not the same amount of rows as Y has observations.", sep = "\n")
  if (is.matrix(vio_space))
    if(any(is.na(vio_space)))
      error_message <- paste(error_message, "There are NA's in vio_space.", sep = "\n")
  if (is.list(vio_space))
    if(any(is.na(unlist(vio_space))))
      error_message <- paste(error_message, "There are NA's in vio_space.", sep = "\n")
  if (alpha > 0.5)
    error_message <- paste(error_message, "alpha cannot be larget than 0.5.", sep = "\n")

  if (!is.null(error_message))
    stop(error_message)

  if (is.null(A1_ind)) A1_ind <- seq_len(NROW(Y))


  Y <- as.matrix(Y)
  D <- as.matrix(D)
  Z <- as.matrix(Z)
  if (!is.null(X)) X <- as.matrix(X)
  n_A1 <- length(A1_ind)

  list_vio_space <- check_vio_space(Z, vio_space)
  vio_space <- list_vio_space$vio_space[A1_ind, ]
  rm_ind <- list_vio_space$rm_ind
  Q <- list_vio_space$Q

  Y_A1 <- Y[A1_ind, ]
  D_A1 <- D[A1_ind, ]
  if (is.null(X)){
    X_A1 <- NULL
  } else {
    X_A1 <- X[A1_ind, ]
  }

  outputs <- tsci_selection(
    Y = Y,
    D = D,
    X = X,
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
  outputs <- append(outputs,
                    list(FirstStage_model = "Specified by User",
                         n_A1 = n_A1,
                         n_A2 = n - n_A1,
                         nsplits = NULL,
                         mult_split_method = NULL,
                         alpha = alpha,
                         ))
  class(outputs) <- c("tsci", "list")
  return(outputs)
}
