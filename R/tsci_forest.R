#' Two Stage Curvature Identification with Random Forests.
#' @description This function implements Two Stage Curvature Identification with the Random Forests. It tests the IV strength and chooses the best violation form, and also constructs the confidence interval for the treatment effect with the selected violation form.
#'
#' @param Y outcome with dimension n by 1.
#' @param D treatment with dimension n by 1.
#' @param Z instrument variable with dimension n by 1.
#' @param X baseline covariates with dimension n by p.
#' @param vio_space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation form selection will be performed according to provided violation forms, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case.
#' @param intercept logic, including the intercept or not in the outcome model, default by TRUE.
#' @param split_prop numeric, proportion of observations used to fit the outcome model.
#' @param num_trees number of trees in Random Forests, default by 200.
#' @param mtry number of covariates to possibly split at in each node of the tree in Random Forests, default by a sequence from round((p+1)/3) to round(2(p+1)/3).
#' @param max_depth maximal tree depth in Random Forests, default by 0, which refers to unlimited depth.
#' @param min_node_size minimal size of each leaf node in Random Forests, default by the set {5, 10, 15}.
#' @param str_thol minimal value of the threshold of IV strength test, default by 10.
#' @param alpha the significance level, default by 0.05.
#' @param parallel One out of \code{"no"}, \code{"multicore"}, or \code{"snow"} specifying the parallelization method used.
#' @param nsplits numeric, number of times the data will be split.
#' @param mult_split_method method to for inference if multi-splitting is performed. Either 'DML' or 'FWER'.
#' @param ncores numeric, the number of cores used. \code{mclapply} form the package \code{parallel} will be called. Parallelization is not supported for Windows.
#' @param cl Either an parallel or snow cluster or \code{NULL}.
#'
#' @return
#' \describe{
#'     \item{\code{Coef_all}}{a series of point estimators of treatment effect corresponding to different violation spaces and the OLS.}
#'     \item{\code{sd_all}}{standard errors of Coef_all.}
#'     \item{\code{CI_all}}{confidence intervals for the treatment effect corresponding to different violation spaces and the OLS.}
#'     \item{\code{Coef_robust}}{the point estimators corresponding to the violation space selected by the robust comparison.}
#'     \item{\code{sd_robust}}{the standard errors of Coef_robust.}
#'     \item{\code{CI_robust}}{confidence intervals for the treatment effect with the violation space selected by the robust comparison.}
#'     \item{\code{iv_str}}{IV strength corresponding to different violation spaces.}
#'     \item{\code{iv_thol}}{the threshold of IV strength test corresponding to different violation spaces.}
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the null violation space and run TSRF only for null violation space. In other cases, violation space selection is performed.}
#'     \item{\code{q_hat}}{the index of estimated violation space corresponding to Qmax.}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid.}
#' }
#' @export
#'
#' @examples
#' # dimension
#' p = 10
#' # sample size
#' n = 1000
#' # interaction value
#' inter.val = 1
#' # the IV strength
#' a = 1
#' # violation strength
#' tau = 1
#' f = function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
#' rho1=0.5
#' # function to generate covariance matrix
#' A1gen=function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'     for(j in 1:p){
#'       A1[i,j]=rho^(abs(i-j))
#'     }
#'   }
#'   A1
#' }
#' Cov=(A1gen(rho1,p+1))
#' mu=rep(0,p+1)
#' # true effect
#' beta=1
#' alpha=as.matrix(rep(-0.3,p))
#' gamma=as.matrix(rep(0.2,p))
#' inter=as.matrix(c(rep(inter.val,5),rep(0,p-5)))
#'
#'
#' # generate the data
#' mu.error=rep(0,2)
#' Cov.error=matrix(c(1,0.5,0.5,1),2,2)
#' Error=MASS::mvrnorm(n, mu.error, Cov.error)
#' W.original=MASS::mvrnorm(n, mu, Cov)
#' W=pnorm(W.original)
#' # instrument variable
#' Z=W[,1]
#' # baseline covariates
#' X=W[,-1]
#' # generate the treatment variable D
#' D=f(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
#' # generate the outcome variable Y
#' Y=D*beta+tau*Z+X%*%gamma+Error[,2]
#'
#'
#' # Two Stage Random Forest
#' output.RF = tsci_forest(Y,D,Z,X)
#' # point estimates
#' output.RF$Coef_robust
#' # standard errors
#' output.RF$sd_robust
#' # confidence intervals
#' output.RF$CI_robust
tsci_forest <- function(Y,
                        D,
                        Z,
                        X,
                        intercept = TRUE,
                        vio_space = NULL,
                        split_prop = 2 / 3,
                        num_trees = NULL,
                        mtry = NULL,
                        max_depth = NULL,
                        min_node_size = NULL,
                        str_thol = 10,
                        alpha = 0.05,
                        parallel = "no",
                        nsplits = 1,
                        mult_split_method = NULL,
                        ncores = 1,
                        cl = NULL) {
  if (!is.null(vio_space)) {
    if (class(vio_space)[1] != "matrix" & class(vio_space)[1] != "list") {
      stop("The violation space must be input as matrix or list")
    }
  }
  if (!is.null(mult_split_method)) {
    if (!(mult_split_method %in% c("FWER", "DML"))) {
      stop("No valid multi-splitting inference method selected. Choose either 'DML' or 'FWER'.")
    }
  }
  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z); X = as.matrix(X);
  n = nrow(X); p = ncol(X)

  do_parallel <- parallelization_setup(parallel = parallel, ncpus = ncores, cl = cl)

  # default value for hyper-parameters
  if (is.null(num_trees)) num_trees <- 200
  if (is.null(mtry)) mtry <- seq(round((p + 1) / 3), round(2 * (p + 1) / 3), by=1)
  if (is.null(max_depth)) max_depth <- 0
  if (is.null(min_node_size)) min_node_size <- c(5, 10, 20)
  if (is.null(nsplits)) nsplits <- 1
  if (is.null(mult_split_method))  mult_split_method <- "DML"

  # grid search
  params_grid <- expand.grid(
    num_trees = num_trees,
    mtry = mtry,
    max_depth = max_depth,
    min_node_size = min_node_size
  )

  # Treatment model fitting
  W <- as.matrix(cbind(Z, X))
  D <- as.matrix(D)
  n <- NROW(W)
  p <- NCOL(W)
  df_treatment <- data.frame(cbind(D, W))
  names(df_treatment) <- c("D", paste("W", seq_len(p), sep = ""))

  # split the data into two parts A1 and A2
  # use A2 to train and use A1 to predict
  n_A1 <- round(split_prop * n)
  n_A2 <- n - n_A1
  A1_ind <- seq_len(n_A1)
  df_treatment_A1 <- df_treatment[A1_ind, ]
  df_treatment_A2 <- df_treatment[-A1_ind, ]

  # Hyperparameter selection
  forest_OOB <- get_forest_parameters(df_treatment_A2 = df_treatment_A2, params_grid = params_grid)

  outputs <- multi_split(df_treatment = df_treatment,
                         Y = Y,
                         D = D,
                         Z = Z,
                         X = X,
                         vio_space = vio_space,
                         A1_ind = A1_ind,
                         intercept = intercept,
                         str_thol = str_thol,
                         alpha = alpha,
                         params = forest_OOB$params_A2,
                         function_hatmatrix = get_forest_hatmatrix,
                         split_prop = split_prop,
                         parallel = parallel,
                         do_parallel = do_parallel,
                         nsplits = nsplits,
                         ncores = ncores,
                         mult_split_method = mult_split_method,
                         cl = cl)

  # Return output
  outputs <- append(outputs, list("mse_cv" = forest_OOB$MSE_CV_A2))
  return(outputs)
}
