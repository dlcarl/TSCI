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
#' @param X observations of baseline covariate(s) used to fit the treatment model. Either a numeric vector of length n
#' or a numeric matrix with dimension n by p or \code{NULL}
#' (if no covariates should be included).
#' @param W (transformed) observations of baseline covariate(s) used to fit the outcome model. Either a numeric vector of length n
#' or a numeric matrix with dimension n by p_w or \code{NULL}
#' (if no covariates should be included).
#' @param vio_space either \code{NULL} or a list with numeric vectors of length n and/or numeric matrices with n rows as elements to
#' specify the violation space candidates.
#' If \code{NULL}, then the violation space candidates are chosen to be a nested sequence
#' of monomials with degree depending on the orders of the polynomials used to fit
#' the treatment model.
#' @param create_nested_sequence logical. If \code{TRUE} the violation space candidates (in form of matrices)
#' are defined sequentially starting with an empty violation matrix and subsequently
#' adding the next element of \code{vio_space} to the current violation matrix.
#' If \code{FALSE} the violation space candidates (in form of matrices) are defined as the elements of \code{vio_space}.
#' See Details for more information.
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
#' @param max_iter number of iterations used in the backfitting algorithm if \code{order_selection_method} is 'backfitting'.
#' @param conv_tol tolerance of convergence in the backfitting algorithm if \code{order_selection_method} is 'backfitting'.
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
#' where \eqn{g(Z_i, X_i)} is estimated using a polynomial basis expansion of the instrumental variables
#' and a linear combination of the baseline covariates and
#' \eqn{h(Z_i X_i)} is approximated using the violation space candidates and by
#' a linear combination the columns in \code{W}. The errors are allowed to be heteroscedastic. \cr \cr
#' The violation space candidates should be in a nested sequence as otherwise nonsensical results can occur.
#' If \code{vio_space} is \code{NULL} the violation space candidates are chosen to be a nested sequence
#' of polynomials of the instrumental variables up to the degrees used to fit the treatment model.
#' This guarantees that the possible spaces the violation lives will be tested.
#' If the functional form of the outcome model is not well-known it is advisable to use the default values
#' for \code{W} and \code{vio_space}. \cr \cr
#'
#' @references
#' \itemize{
#' \item{Zijian Guo, and Peter Bühlmann. Two Stage Curvature Identification with
#' Machine Learning: Causal Inference with Possibly Invalid Instrumental Variables.
#' \emph{arXiv:2203.12808}, 2022}
#' }
#'
#' @seealso
#' \code{\link[TSML]{tsci_forest}} for TSCI with random forest. \cr \cr
#' \code{\link[TSML]{tsci_boosting}} for TSCI with boosting. \cr \cr
#' \code{\link[TSML]{tsci_secondstage}} for TSCI with user provided hat matrix. \cr \cr
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
                      W = X,
                      vio_space = NULL,
                      create_nested_sequence = TRUE,
                      intercept = TRUE,
                      min_order = 1,
                      max_order = 10,
                      exact_order = NULL,
                      order_selection_method = c("backfitting", "grid search"),
                      max_iter = 100,
                      conv_tol = 10^-6,
                      gcv = T,
                      nfolds = 5,
                      str_thol = 10,
                      alpha = 0.05) {

  # checks that input is in the correct format
  check_input(Y = Y,
              D = D,
              Z = Z,
              X = X,
              W = W,
              vio_space = vio_space,
              create_nested_sequence = create_nested_sequence,
              intercept = intercept,
              min_order = min_order,
              max_order = max_order,
              exact_order = exact_order,
              max_iter = max_iter,
              conv_tol = conv_tol,
              gcv = gcv,
              nfolds = nfolds,
              str_thol = str_thol,
              alpha = alpha,
              tsci_method = "poly")
  order_selection_method <- match.arg(order_selection_method)


  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z)
  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(W)) W <- as.matrix(W)

  n <- NROW(Y)
  p <- NCOL(Z) + ifelse(is.null(X), 0, NCOL(X))

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
  df_treatment <- data.frame(cbind(D, Z, X))
  names(df_treatment) <- c("D", paste("B", seq_len(p), sep = ""))

  # hyperparameter tuning
  poly_CV <- get_poly_parameters(df_treatment = df_treatment,
                                 params_list = params_list,
                                 order_selection_method = order_selection_method,
                                 max_iter = max_iter,
                                 conv_tol,
                                 gcv = gcv,
                                 nfolds = nfolds)

  if (is.null(vio_space)) {
    vio_space <- create_monomials(Z = Z,
                                  degree = unlist(poly_CV$params)[seq_len(NCOL(Z))],
                                  type = "monomials_main")
    create_nested_sequence <- TRUE

  }
  list_vio_space <- build_vio_space_candidates(vio_space = vio_space,
                                               create_nested_sequence = create_nested_sequence)

  if (!(list_vio_space$nested_sequence))
    warning("Sequence of violation space candidates is not nested. Results might be nonsensical.")

  outputs <- tsci_fit(df_treatment = df_treatment,
                      Y = Y,
                      D = D,
                      Z = Z,
                      W = W,
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
