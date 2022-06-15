#' TSCI with multiple data splitting.
#'
#' @param df_treatment a data frame containing the subset of observations of the treatment model.
#' @param Y outcome vector with dimension n by 1.
#' @param D treatment vector with dimension n by 1.
#' @param Z instrument variable with dimension n by 1.
#' @param X baseline covariates with dimension n by p.
#' @param vio_space vio_space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation form selection will be performed according to provided violation forms, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case.
#' @param A1_ind the indices of samples in A1 in the first split.
#' @param intercept logical, including the intercept or not in the outcome model, default by TRUE.
#' @param str_thol minimal value of the threshold of IV strength test, default by 10.
#' @param alpha the significance level, default by 0.05.
#' @param params a list containing the hyperparameters of the treatment model fitting method.
#' @param function_hatmatrix a function to get the hat matrix of the treatment model.
#' @param split_prop numeric, proportion of observations used to fit the outcome model.
#' @param nsplits numeric, number of times the data will be split.
#' @param parallel character, one out of \code{"no"}, \code{"multicore"}, or \code{"snow"} specifying the parallelization method used.
#' @param do_parallel logical, specifies if parallel computing should be performed.
#' @param ncores numeric, the number of cores used if multi_splitting is \code{TRUE}. \code{mclapply} form the package \code{parallel} will be called. Parallelization is not supported for Windows.
#' @param mult_split_method method to for inference if multi-splitting is performed. Either 'DML' or 'FWER'.
#' @param cl Either an parallel or snow cluster or \code{NULL}.
#'
#' @return
#'     \item{\code{Coef_all}}{the median over the multiple data splits of a series of point estimators of treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{sd_all}}{standard errors of Coef_all}
#'     \item{\code{CI_all}}{confidence intervals for the treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{Coef_robust}}{the median over the multiple data splits of the point estimators corresponding to the violation space selected by the robust comparison}
#'     \item{\code{sd_robust}}{the standard errors of Coef_robust}
#'     \item{\code{CI_robust}}{confidence intervals for the treatment effect with the violation space selected by the robust comparison}
#'     \item{\code{iv_str}}{the median over the multiple data splits of IV strength corresponding to different violation spaces}
#'     \item{\code{iv_thol}}{the median over the multiple data splits of the threshold of IV strength test corresponding to different violation spaces}
#'     \item{\code{Qmax}}{the median over the multiple data splits of the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the first violation space and run TSRF only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q_hat}}{the median over the multiple data splits of the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{the median over the multiple data splits of invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' @noRd
#'
multi_split <- function(df_treatment,
                        Y,
                        D,
                        Z,
                        X,
                        vio_space,
                        A1_ind,
                        intercept,
                        str_thol,
                        alpha,
                        params,
                        function_hatmatrix,
                        split_prop,
                        parallel,
                        do_parallel,
                        nsplits,
                        ncores,
                        mult_split_method,
                        cl) {

  if (is.null(vio_space)) {
    Q <- 4
  }

  if (!is.null(vio_space)) {
    if (class(vio_space)[1] == "list") {
      Q <- length(vio_space) + 1
    } else if (class(vio_space)[1] == "matrix") {
      Q <- ncol(vio_space) + 1
    }
  }

  tsci_parallel <- local({
    df_treatment
    Y
    D
    Z
    X
    vio_space
    A1_ind
    intercept
    str_thol
    alpha
    params
    function_hatmatrix
    ncores
    Q
    function(colnames.cluster) {
      tryCatch_WEM(tsci_fit(
        df_treatment = df_treatment,
        Y = Y,
        D = D,
        Z = Z,
        X = X,
        vio_space = vio_space,
        intercept = intercept,
        str_thol = str_thol,
        split_prop = split_prop,
        alpha = alpha,
        params = params,
        function_hatmatrix = function_hatmatrix
      ), tsci_fit_NA_return(Q = Q))}
  })

  if (do_parallel) {
    if (parallel == "multicore") {
      list_outputs <- parallel::mclapply(seq_len(nsplits), tsci_parallel, mc.cores = ncores)
    } else if (parallel == "snow") {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncores))
        # export the namespace of dmlalg in order for the use the functions
        # of the package TSML on the workers
        parallel::clusterExport(cl, varlist = getNamespaceExports("TSML"))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        list_outputs <- parallel::parLapply(cl, seq_len(nsplits), tsci_parallel)
        parallel::stopCluster(cl)
        cl <- NULL # overwrite object which is responsible for the connection
      } else list_outputs <- parallel::parLapply(cl, seq_len(nsplits), tsci_parallel)
    }
  } else list_outputs <- lapply(seq_len(nsplits), tsci_parallel)

  Coef_all_matrix <-
    matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$Coef_all), use.names = FALSE),
           ncol = length(list_outputs[[1]]$value$Coef_all), byrow = TRUE)
  colnames(Coef_all_matrix) <- names(list_outputs[[1]]$value$Coef_all)
  sd_all_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$sd_all),
                                 use.names = FALSE), ncol = length(list_outputs[[1]]$value$sd_all), byrow = TRUE)
  colnames(sd_all_matrix) <- names(list_outputs[[1]]$value$sd_all)
  Coef_robust_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$Coef_robust),
                                      use.names = FALSE), ncol = length(list_outputs[[1]]$value$Coef_robust), byrow = TRUE)
  colnames(Coef_robust_matrix) <- names(list_outputs[[1]]$value$Coef_robust)
  sd_robust_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$sd_robust),
                                    use.names = FALSE), ncol = length(list_outputs[[1]]$value$sd_robust), byrow = TRUE)
  colnames(sd_robust_matrix) <- names(list_outputs[[1]]$value$sd_robust)
  iv_str_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$iv_str),
                                 use.names = FALSE), ncol = length(list_outputs[[1]]$value$iv_str), byrow = TRUE)
  colnames(iv_str_matrix) <- names(list_outputs[[1]]$value$iv_str)
  iv_thol_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$iv_thol),
                                  use.names = FALSE), ncol = length(list_outputs[[1]]$value$iv_thol), byrow = TRUE)
  colnames(iv_thol_matrix) <- names(list_outputs[[1]]$value$iv_thol)
  Qmax_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$Qmax),
                               use.names = FALSE), ncol = length(list_outputs[[1]]$value$Qmax), byrow = TRUE)
  colnames(Qmax_matrix) <- names(list_outputs[[1]]$value$Qmax)
  q_hat_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$q_hat),
                                use.names = FALSE), ncol = length(list_outputs[[1]]$value$q_hat), byrow = TRUE)
  colnames(q_hat_matrix) <- names(list_outputs[[1]]$value$q_hat)
  invalidity_matrix <- matrix(unlist(lapply(list_outputs, FUN = function(x) x$value$invalidity),
                                     use.names = FALSE), ncol = length(list_outputs[[1]]$value$invalidity), byrow = TRUE)
  colnames(invalidity_matrix) <- names(list_outputs[[1]]$value$invalidity)

  Coef_all <- apply(Coef_all_matrix, 2, FUN = stats::median)
  Coef_robust <- apply(Coef_robust_matrix, 2, FUN = stats::median)
  iv_str <- apply(iv_str_matrix, 2, FUN = stats::median)
  iv_thol <- apply(iv_thol_matrix, 2, FUN = stats::median)
  Qmax <- apply(Qmax_matrix, 2, FUN = stats::median)
  q_hat <- apply(q_hat_matrix, 2, FUN = stats::median)
  invalidity <- apply(invalidity_matrix, 2, FUN = stats::median)

  if (mult_split_method == "FWER") {
    CI_all <- sapply(seq_len(NCOL(Coef_all_matrix)),
      FUN = function(j) {
        beta <- Coef_all_matrix[, j]
        se <- sd_all_matrix[, j]
        lower <- stats::median(beta - 1 * stats::qnorm(1 - alpha / 4) * se)
        upper <- stats::median(beta + 1 * stats::qnorm(1 - alpha / 4) * se)
        return(c(lower, upper))
      }
    )
    rownames(CI_all) <- c("lower", "upper")
    colnames(CI_all) <- colnames(Coef_all_matrix)

    CI_robust <- sapply(seq_len(NCOL(Coef_robust_matrix)),
      FUN = function(j) {
        beta <- Coef_robust_matrix[, j]
        se <- sd_robust_matrix[, j]
        lower <- stats::median(beta - 1 * stats::qnorm(1 - alpha / 4) * se)
        upper <- stats::median(beta + 1 * stats::qnorm(1 - alpha / 4) * se)
        return(c(lower, upper))
      }
    )
    rownames(CI_robust) <- c("lower", "upper")
    colnames(CI_robust) <- colnames(Coef_robust_matrix)

    # this method does not provide standard error estimates.
    sd_all <- rep("NA", times = NCOL(sd_all_matrix))
    names(sd_all) <- colnames(sd_all_matrix)
    sd_robust <- rep("NA", times = NCOL(sd_robust_matrix))
    names(sd_robust) <- colnames(sd_robust_matrix)
  } else if (mult_split_method == "DML") {
    sd_all <- sapply(seq_len(NCOL(sd_all_matrix)),
      FUN = function(j) {
        stats::median(sqrt(sd_all_matrix[, j]^2 +
          (Coef_all_matrix[, j] - stats::median(Coef_all_matrix[, j]))^2))
      }
    )
    names(sd_all) <- colnames(sd_all_matrix)
    CI_all <- sapply(seq_len(NCOL(Coef_all_matrix)),
      FUN = function(j) {
        lower <- stats::median(Coef_all_matrix[, j]) -
          stats::qnorm(1 - alpha / 2) * sd_all[j]
        upper <- stats::median(Coef_all_matrix[, j]) +
          stats::qnorm(1 - alpha / 2) * sd_all[j]
        return(c(lower, upper))
      }
    )
    rownames(CI_all) <- c("lower", "upper")
    colnames(CI_all) <- colnames(Coef_all_matrix)

    sd_robust <- sapply(seq_len(NCOL(sd_robust_matrix)),
      FUN = function(j) {
        stats::median(sqrt(sd_robust_matrix[, j]^2 +
          (Coef_robust_matrix[, j] - stats::median(Coef_robust_matrix[, j]))^2))
      }
    )
    names(sd_robust) <- colnames(sd_robust_matrix)
    CI_robust <- sapply(seq_len(NCOL(Coef_robust_matrix)),
      FUN = function(j) {
        lower <- stats::median(Coef_robust_matrix[, j]) -
          stats::qnorm(1 - alpha / 2) * sd_robust[j]
        upper <- stats::median(Coef_robust_matrix[, j]) +
          stats::qnorm(1 - alpha / 2) * sd_robust[j]
        return(c(lower, upper))
      }
    )
    rownames(CI_robust) <- c("lower", "upper")
    colnames(CI_robust) <- colnames(Coef_robust_matrix)
  } else {
    # if no valid inference method was selected, remove measures of uncertainty.
    sd_all <- rep("NA", times = NCOL(sd_all_matrix))
    names(sd_all) <- colnames(sd_all_matrix)
    CI_all <- matrix("NA", nrow = 2, ncol = NCOL(Coef_all_matrix))
    rownames(CI_all) <- c("lower", "upper")
    colnames(CI_all) <- colnames(Coef_all_matrix)
    sd_robust <- rep("NA", times = NCOL(sd_robust_matrix))
    names(sd_robust) <- colnames(sd_robust_matrix)
    CI_robust <- matrix("NA", nrow = 2, ncol = NCOL(Coef_robust_matrix))
    rownames(CI_robust) <- c("lower", "upper")
    colnames(CI_robust) <- colnames(Coef_all_matrix)
  }
  returnList <- list(
    Coef_all = Coef_all,
    sd_all = sd_all,
    CI_all = CI_all,
    Coef_robust = Coef_robust,
    sd_robust = sd_robust,
    CI_robust = CI_robust,
    iv_str = iv_str,
    iv_thol = iv_thol,
    Qmax = Qmax,
    q_hat = q_hat,
    invalidity = invalidity
  )
  returnList
}
