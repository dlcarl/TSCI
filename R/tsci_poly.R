#' Two Stage Curvature Identification with Polynomial Basis Expansion
#' @description \code{tsci_poly} implements Two Stage Curvature Identification
#' (Guo and Bühlmann 2022) with a basis expansion by monomials. Through a data-dependent way it
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
#' If \code{NULL}, then the violation space candidates are chosen to be a nested sequence
#' of monomials with degree depending on the orders of the polynomials used to fit
#' the treatment model.
#' @param intercept logical. If \code{TRUE} an intercept is included in the outcome model.
#' @param min_order either a single numeric value or a numeric vector of length s specifying
#' the smallest order of polynomials to use in the selection of the treatment model. If a
#' single numeric value, all the polynomials of all instrumental variables use this value.
#' @param max_order either a single numeric value or a numeric vector of length s specifying
#' the largest order of polynomials to use in the selection of the treatment model. If a
#' single numeric value, all the polynomials of all instrumental variables use this value.
#' @param exact_order either a single numeric value or a numeric vector of length s specifying
#' the exact order of polynomials to use in the treatment model. If a
#' single numeric value, all the polynomials of all instrumental variables use this value.
#' @param order_selection_method method used to select the best fitting order of polynomials
#' for the treatment model. Must be either 'backfitting' or 'grid search'.
#' 'grid search' can be very slow if the number of instruments is large.
#' @param gcv logical. If \code{TRUE} the generalized cross validation mean squared error is used
#' to determine the best fitting order of polynomials for the treatment model.
#' If \code{FALSE} k-fold cross validation is used instead.
#' @param nfolds number of folds used for the k-fold cross-validation if \code{gcv} is \code{FALSE}.
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
#' where \eqn{g(Z_i, X_i)} is estimated using L2 boosting with regression trees as base learners and
#' \eqn{h(Z_i X_i)} is approximated using the violation space candidates and by
#' a linear combination of baseline covariates. The errors are allowed to be heteroscedastic.
#' The violation space candidates are required to be in a nested sequence. The specification
#' of suitable violation space candidates is a crucial step because a poor approximation
#' of \eqn{h(Z_i, X_i)} might not address the bias caused by the violation of the IV assumption sufficiently.
#' If \code{vio_space} is \code{NULL} the violation space candidates are chosen to be
#' \eqn{{}, {Z1, Z2, ..., Zs}, {Z1, Z2, ..., Zs, Z1^2, Z2^2, ..., Z2^2}} and
#' \eqn{{Z1, Z2, ..., Zs, Z1^2, Z2^2, ..., Z2^2, Z1^3, Z2^3, ..., Z2^3}} thus implicitly assuming
#' that there are no interactions between the instruments and the covariates and
#' between the instruments themselves in the outcome model. \cr \cr
#'
#' @references
#' \itemize{
#' \item{Zijian Guo, and Peter Bühlmann. Two Stage Curvature Identification with
#' Machine Learning: Causal Inference with Possibly Invalid Instrumental Variables.
#' \emph{arXiv:2203.12808}, 2022}
#' }
#'
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
#'
#' # Two Stage Polynomials
#' output_PO <- tsci_poly(Y, D, Z, X)
#' # point estimates
#' output_PO$Coef_robust
#' # standard errors
#' output_PO$sd_robust
#' # confidence intervals
#' output_PO$CI_robust
tsci_poly <- function(Y,
                      D,
                      Z,
                      X = NULL,
                      vio_space = NULL,
                      intercept = TRUE,
                      min_order = 1,
                      max_order = 10,
                      exact_order = NULL,
                      order_selection_method = "backfitting",
                      gcv = T,
                      nfolds = 5,
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
    error_message <- paste(error_message, "X is not numeric.", sep = "\n")
  if (!is.numeric(min_order) & !is.null(min_order))
    error_message <- paste(error_message, "min_order is not numeric.", sep = "\n")
  if (!is.numeric(max_order) & !is.null(max_order))
    error_message <- paste(error_message, "max_order is not numeric.", sep = "\n")
  if (!is.numeric(exact_order) & !is.null(exact_order))
    error_message <- paste(error_message, "exact_order is not numeric.", sep = "\n")
  if (!is.logical(gcv))
    error_message <- paste(error_message, "gcv is neither TRUE nor FALSE.", sep = "\n")
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
  if (is.matrix(vio_space))
    if(any(is.na(vio_space)))
      error_message <- paste(error_message, "There are NA's in vio_space.", sep = "\n")
  if (is.list(vio_space))
    if(any(is.na(unlist(vio_space))))
      error_message <- paste(error_message, "There are NA's in vio_space.", sep = "\n")
  if (alpha > 1)
    error_message <- paste(error_message, "alpha cannot be larger than 1.", sep = "\n")
  if (is.null(exact_order) & (is.null(min_order) | is.null(max_order)))
    error_message <- paste(error_message, "Either exact_order or min_order and max_order must be specified.", sep = "\n")
  if (length(max_order) != length(min_order))
    error_message <- paste(error_message, "min_order and max_order must be of same length.", sep = "\n")
  if (length(max_order) > 1 & length(max_order) != NCOL(Z))
    error_message <- paste(error_message, "min_order and max_order have invalid length.", sep = "\n")
  if (!is.null(min_order) & any(min_order < 1))
    error_message <- paste(error_message, "Values in min_order cannot be smaller than 1.", sep = "\n")
  if (!is.null(max_order) & any(max_order < min_order))
    error_message <- paste(error_message, "Values in max_order cannot be smaller than values in min_order.", sep = "\n")
  if (!(order_selection_method %in% c("backfitting", "grid search")))
    error_message <- paste(error_message, "No valid order selection method
                           selected. Choose either 'backfitting' or 'grid search'.", sep = "\n")
  if (!is.null(error_message))
    stop(error_message)


  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z)
  if (!is.null(X)) X <- as.matrix(X)


  # grid search
  if (is.null(exact_order)) {
    if (length(max_order) == 1) max_order <- rep(max_order, NCOL(Z))
    if (length(min_order) == 1) min_order <- rep(min_order, NCOL(Z))
  } else if (!is.null(exact_order)) {
    if (length(exact_order) == 1) {
      max_order <- rep(exact_order, NCOL(Z))
      min_order <- rep(exact_order, NCOL(Z))
    } else {
      max_order <- exact_order
      min_order <- exact_order
    }
  }

  if (!is.null(X)) {
    max_order <- c(max_order, rep(1, NCOL(X)))
    min_order <- c(min_order, rep(1, NCOL(X)))
  }

  params_list <- lapply(seq_len(p), FUN = function(i) seq(min_order[i], max_order[i], by = 1))

  # Treatment model fitting
  W <- as.matrix(cbind(Z, X))
  df_treatment <- data.frame(cbind(D, W))
  names(df_treatment) <- c("D", paste("W", seq_len(p), sep = ""))

  # hyperparameter tuning
  poly_CV <- get_poly_parameters(df_treatment = df_treatment,
                                 params_list = params_list,
                                 order_selection_method = order_selection_method,
                                 gcv = gcv,
                                 nfolds = nfolds)

  if (is.null(vio_space)) vio_space <- create_monomials(Z = Z,
                                                        degree = max(unlist(poly_CV$params)[seq_len(NCOL(Z))]),
                                                        type = "monomials_main")
  list_vio_space <- build_vio_space_candidates(Z, vio_space)

  outputs <- tsci_fit(df_treatment = df_treatment,
                      Y = Y,
                      D = D,
                      Z = Z,
                      X = X,
                      list_vio_space = list_vio_space,
                      intercept = intercept,
                      str_thol = str_thol,
                      split_prop = 1,
                      alpha = alpha,
                      params = poly_CV$params,
                      function_hatmatrix = get_poly_hatmatrix)

  # Return output
  outputs <- append(outputs,
                    list(mse = poly_CV$mse,
                         FirstStage_model = "OLS with Polynomials",
                         FirstStage_params = poly_CV$params,
                         n_A1 = n,
                         n_A2 = 0,
                         nsplits = 0,
                         mult_split_method = "No sample splitting was performed",
                         alpha = alpha))
  class(outputs) <- c("tsci", "list")
  return(outputs)
}
