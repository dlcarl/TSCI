#' Two Stage Curvature Identification with Basis Splines
#' @description \code{tsci_boosting} implements Two Stage Curvature Identification
#' (Guo and Bühlmann 2022) with basis splines. Through a data-dependent way it
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
#' If \code{NULL}, then the violation space candidates are the spaces of polynomials up to the 3-th order.
#' See Details for more information.
#' @param intercept logical. If \code{TRUE} an intercept is included in the outcome model.
#' @param norder xxx
#' @param nfolds number of folds used for cross-validation to choose best parameter combination.
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
#' To avoid overfitting bias the data is randomly split into two subsets \eqn{A1} and \eqn{A2}
#' where the proportion of number of observations in the two sets is specified by \code{split_prop}.
#' \eqn{A2} is used to train the boosting model and \eqn{A1} is used to fit the outcome model. \cr \cr
#' The package \code{xgboost} is used for boosting. If any of \code{nrounds},
#' \code{eta}, \code{max_depth}, \code{subsample} or \code{colsample_bytree} has more than one value,
#' the best parameter combination is chosen by minimizing the cross-validation mean squared error. \cr \cr
#' The violation space candidates are required to be in a nested sequence. The specification
#' of suitable violation space candidates is a crucial step because a poor approximation
#' of \eqn{h(Z_i, X_i)} might not address the bias caused by the violation of the IV assumption sufficiently.
#' If \code{vio_space} is \code{NULL} the violation space candidates are chosen to be
#' \eqn{{}, {Z1, Z2, ..., Zs}, {Z1, Z2, ..., Zs, Z1^2, Z2^2, ..., Z2^2}} and
#' \eqn{{Z1, Z2, ..., Zs, Z1^2, Z2^2, ..., Z2^2, Z1^3, Z2^3, ..., Z2^3}} thus implicitly assuming
#' that there are no interactions between the instruments and the covariates and
#' between the instruments themselves in the outcome model. \cr \cr
#' If \code{nsplits} is larger than 1, point estimates are aggregated by medians
#' and standard errors, p-values and confidence intervals are obtained by the method
#' specified by the parameter \code{mult_split_method}. 'DML' uses the approach by
#' Chernozhukov et al. (2018). 'FWER' uses the approach by Meinshausen et al. (2009)
#' and controls for the family-wise error rate. 'FWER' does not provide standard errors.
#' For large sample sizes a large values for \code{nsplits} can lead to a high
#' running time as for each split a new hat matrix must be calculated.
#' The same parameter combination for the boosting model is used for all n splits and
#' is selected by a separate data split.
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
                         norder = seq(1, 6, by = 1),
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
    error_message <- paste(error_message, "alpha cannot be larget than 1.", sep = "\n")
  if (!is.null(error_message))
    stop(error_message)


  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z)
  if (!is.null(X)) X <- as.matrix(X)


  # grid search
  params_grid <- expand.grid(
    norder = norder,
    placeholder = 1
  )

  # Treatment model fitting
  W <- as.matrix(cbind(Z, X))
  df_treatment <- data.frame(cbind(D, W))
  names(df_treatment) <- c("D", paste("W", seq_len(p), sep = ""))

  # hyperparameter tuning
  poly_CV <- get_poly_parameters(df_treatment, params_grid = params_grid, nfolds = nfolds)

  if (is.null(vio_space)) vio_space <- poly(as.vector(Z), degree = poly_CV$params$norder, raw = TRUE, simple = TRUE)
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
