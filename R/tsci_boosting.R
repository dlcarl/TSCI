#' Two Stage Curvature Identification with Boosting
#' @description \code{tsci_boosting} implements Two Stage Curvature Identification
#' (Guo and Bühlmann 2022) with boosting. Through a data-dependent way it
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
#' @param vio_space  list with numeric vectors of length n and/or numeric matrices with n rows as elements to
#' specify the violation space candidates. See Details for more information.
#' @param create_nested_sequence logical. If \code{TRUE} the violation space candidates (in form of matrices)
#' are defined sequentially starting with an empty violation matrix and subsequently
#' adding the next element of \code{vio_space} to the current violation matrix.
#' If \code{FALSE} the violation space candidates (in form of matrices) are defined as the elements of \code{vio_space}.
#' See Details for more information.
#' @param intercept logical. If \code{TRUE} an intercept is included in the outcome model.
#' @param split_prop proportion of observations used to fit the outcome model. Has to be a value in (0, 1).
#' @param nrounds number of boosting iterations.
#' Can either be a single integer value or a vector containing multiple integer values to try.
#' @param eta learning rate of the boosting algorithm.
#' Can either be a single numeric value or a vector containing multiple numeric values to try.
#' @param max_depth maximal tree depth.
#' Can either be a single integer value or a vector containing multiple integer values to try.
#' @param subsample subsample ratio of the training instance.
#' Can either be a single numeric value or a vector containing multiple numeric values to try.
#' Has to be a value in (0, 1].
#' @param colsample_bytree subsample ratio of columns when constructing each tree.
#' Can either be a single numeric value or a vector containing multiple numeric values to try.
#' Has to be a value in (0, 1].
#' @param early_stopping logical. If \code{TRUE} early stopping will be applied to
#' choose the optimal number of boosting iteration using the cross-validation mean squared error.
#' @param nfolds number of folds used for cross-validation to choose best parameter combination.
#' @param str_thol minimal value of the threshold of IV strength test.
#' @param alpha the significance level.
#' @param nsplits number of times the data will be split. Has to be an integer larger or equal 1.
#' @param mult_split_method method to calculate the standard errors and p-values.
#' and to construct the confidence intervals if multi-splitting is performed. Default is "FWER" if \code{nsplits} > 1 and
#' "DML" otherwise. See Details.
#' @param parallel One out of \code{"no"}, \code{"multicore"}, or \code{"snow"} specifying the parallelization method used.
#' @param ncores the number of cores to use.
#' @param cl either an parallel or snow cluster or \code{NULL}.
#' @param raw_output logical. If \code{TRUE} the coefficient and standard error estimates of each split will be returned.
#' This is only needed if \code{mult_split_method} equals "FWER" and the function \code{confint} will be used. Default is
#' \code{TRUE} if \code{mult_split_method} is \code{TRUE} and \code{FALSE} otherwise.
#' @param B number of bootstrap samples.
#' Bootstrap methods are used to calculate the iv strength threshold and in the violation space selection.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'     \item{\code{Coef_all}}{a series of point estimates of the treatment effect
#'     for the different violation space candidates and the OLS estimate.}
#'     \item{\code{sd_all}}{standard errors of Coef_all.}
#'     \item{\code{pval_all}}{p-values of the treatment effect estimates for the
#'     different violation space candidates and for the OLS estimate.}
#'     \item{\code{CI_all}}{confidence intervals for the treatment effect for the
#'     different violation space candidates and for the OLS estimate.}
#'     \item{\code{Coef_robust}}{the point estimators of the treatment effect for
#'     the selected violation spaces.}
#'     \item{\code{sd_robust}}{the standard errors of Coef_robust.}
#'     \item{\code{pval_robust}}{p-values of the treatment effect estimates for the
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
#'     \item{\code{FirstStage_rse}}{residual standard error of the fitted treatment model.}
#'     \item{\code{FirstStage_Rsquared}}{R-squared of the fitted treatment model.}
#'     \item{\code{FirstStage_rse}}{a named vector containing the residual standard error of the fitted outcome model for each violation space candidate.}
#'     \item{\code{FirstStage_Rsquared}}{a named vector containing the R-squared of the fitted outcome model for each violation space candidate.}
#'     \item{\code{mse}}{the out-of-sample mean squared error of the fitted treatment model.}
#'     \item{\code{FirstStage_model}}{the method used to fit the treatment model.}
#'     \item{\code{FirstStage_params}}{the hyperparameter combination used to fit the treatment model.}
#'     \item{\code{n_A1}}{number of observations in A1.}
#'     \item{\code{n_A2}}{number of observations in A2.}
#'     \item{\code{nsplits}}{number of data splits performed.}
#'     \item{\code{mult_split_method}}{the method used to calculate the standard errors and p-values.}
#'     \item{\code{alpha}}{the significance level used.}
#'}
#'
#' @details The treatment and outcome models are assumed to be of the following forms:
#' \deqn{D_i = f(Z_i, X_i) + \delta_i}
#' \deqn{Y_i = \beta * D_i + g(Z_i, X_i) + \epsilon_i}
#' where \eqn{f(Z_i, X_i)} is estimated using L2 boosting with regression trees as base learners and
#' \eqn{g(Z_i X_i)} is approximated using the violation space candidates and by
#' a linear combination of the columns in \code{W}. The errors are allowed to be heteroscedastic.
#' To avoid overfitting bias the data is randomly split into two subsets \eqn{A1} and \eqn{A2}
#' where the proportion of number of observations in the two sets is specified by \code{split_prop}.
#' \eqn{A2} is used to train the boosting model and \eqn{A1} is used to fit the outcome model. \cr \cr
#' The package \code{xgboost} is used for boosting. If any of \code{nrounds},
#' \code{eta}, \code{max_depth}, \code{subsample} or \code{colsample_bytree} has more than one value,
#' the best parameter combination is chosen by minimizing the cross-validation mean squared error. \cr \cr
#' The violation space candidates should be in a nested sequence as otherwise nonsensical results can occur. The specification
#' of suitable violation space candidates is a crucial step because a poor approximation
#' of \eqn{g(Z_i, X_i)} might not address the bias caused by the violation of the IV assumption sufficiently.
#' The function \code{\link[TSML]{create_monomials}} can be used to create a predefined sequence of violation space candidates (monomials).  \cr \cr
#' \code{W} should be chosen to be flexible enough to approximate the functional form of how the covariates affect the outcome well
#' as otherwise the treatment estimator might be biased. \cr \cr
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
#' @seealso
#' \code{\link[TSML]{tsci_forest}} for TSCI with random forest. \cr \cr
#' \code{\link[TSML]{tsci_poly}} for TSCI with polynomial basis expansion. \cr \cr
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
#' # Two Stage L2 Boosting
#' vio_space <- create_monomials(Z, 4, "monomials_main")
#' output_BO <- tsci_boosting(Y, D, Z, X, vio_space = vio_space)
#' # point estimates
#' output_BO$Coef_robust
#' # standard errors
#' output_BO$sd_robust
#' # confidence intervals
#' output_BO$CI_robust
tsci_boosting <- function(Y,
                          D,
                          Z,
                          X = NULL,
                          W = X,
                          vio_space,
                          create_nested_sequence = TRUE,
                          intercept = TRUE,
                          split_prop = 2 / 3,
                          nrounds = 50,
                          eta = 0.3,
                          max_depth = 6,
                          subsample = 1,
                          colsample_bytree = 1,
                          early_stopping = TRUE,
                          nfolds = 5,
                          str_thol = 10,
                          alpha = 0.05,
                          parallel = "no",
                          nsplits = 10,
                          mult_split_method = c("FWER", "DML"),
                          ncores = 1,
                          cl = NULL,
                          raw_output = NULL,
                          B = 300) {

  # checks that input is in the correct format
  check_input(Y = Y,
              D = D,
              Z = Z,
              X = X,
              W = W,
              vio_space = vio_space,
              create_nested_sequence = create_nested_sequence,
              intercept = intercept,
              split_prop = split_prop,
              nrounds = nrounds,
              eta = eta,
              max_depth = max_depth,
              subsample = subsample,
              colsample_bytree = colsample_bytree,
              early_stopping = early_stopping,
              nfolds = nfolds,
              str_thol = str_thol,
              alpha = alpha,
              parallel = parallel,
              nsplits = nsplits,
              ncores = ncores,
              cl = cl,
              raw_output = raw_output,
              B = B,
              tsci_method = "boosting")
  # if data is split only once then the default aggregation method is "DML" as
  # this gives the same results as not aggregating. Otherwise "FWER" is used
  # because it controls for the family-wise error rate
  if (missing(mult_split_method)) {
    mult_split_method <- ifelse(nsplits > 1, "FWER", "DML")
  }
  mult_split_method <- match.arg(mult_split_method)
  # if TRUE returns the estimate of the treatment effect and its standard error
  # for each data split. Is needed for the confint method if mult_split_method is  "FWER"
  if (is.null(raw_output)) {
    raw_output <- ifelse(mult_split_method == "FWER", TRUE, FALSE)
  }

  # stores variables as matrices as matrix multiplications will be performed later
  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z)
  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(W)) W <- as.matrix(W)

  n <- NROW(Y)
  p <- NCOL(Z) + ifelse(is.null(X), 0, NCOL(X))

  # initializes parallelization setup
  do_parallel <- parallelization_setup(parallel = parallel, ncpus = ncores, cl = cl)

  # sets up grid search over the hyperparameter combinations
  params_grid <- expand.grid(
    nrounds = nrounds,
    eta = eta,
    max_depth = max_depth,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    early_stopping = early_stopping,
    lambda = 0
  )

  # creates the dataframe used to fit the treatment model
  df_treatment <- data.frame(cbind(D, Z, X))
  names(df_treatment) <- c("D", paste("B", seq_len(p), sep = ""))

  # splits the data into two parts A1 and A2
  # A2 will be used to train the treatment model and the hat matrix will be calculated for A1
  n_A1 <- round(split_prop * n)
  n_A2 <- n - n_A1
  A1_ind <- sample(seq_len(n), n_A1)
  df_treatment_A1 <- df_treatment[A1_ind, ]
  df_treatment_A2 <- df_treatment[-A1_ind, ]

  # grid search through params_grid to identify the hyperparameter combination that minimizes the CV MSE
  treeboost_CV <- get_l2boost_parameters(
    df_treatment_A2 = df_treatment_A2,
    params_grid = params_grid,
    nfolds = nfolds)

  # calls tsci_multisplit which splits the data n_splits time into A1 and A2.
  # fits the treatment model with A2, calculates the hat matrix for A1
  # and subsequently performs violation space selection and treatment effect estimation.
  outputs <- tsci_multisplit(df_treatment = df_treatment,
                             Y = Y,
                             D = D,
                             Z = Z,
                             W = W,
                             vio_space = vio_space,
                             create_nested_sequence = create_nested_sequence,
                             A1_ind = A1_ind,
                             intercept = intercept,
                             str_thol = str_thol,
                             alpha = alpha,
                             params = treeboost_CV$params_A2,
                             function_hatmatrix = get_l2boost_hatmatrix,
                             split_prop = split_prop,
                             parallel = parallel,
                             do_parallel = do_parallel,
                             nsplits = nsplits,
                             ncores = ncores,
                             mult_split_method = mult_split_method,
                             cl = cl,
                             raw_output = raw_output,
                             B = B)

  # returns output
  outputs <- append(outputs,
                    list(mse = treeboost_CV$mse,
                         FirstStage_model = "L2 Gradient Tree Boosting",
                         FirstStage_params = treeboost_CV$params_A2,
                         n_A1 = n_A1,
                         n_A2 = n_A2,
                         nsplits = nsplits,
                         mult_split_method = mult_split_method,
                         alpha = alpha))
  class(outputs) <- c("tsci", "list")
  return(outputs)
}
