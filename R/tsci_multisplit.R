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
tsci_multisplit <- function(df_treatment,
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
  list_vio_space <- check_vio_space(Z, vio_space)

  tsci_parallel <- local({
    df_treatment
    Y
    D
    Z
    X
    list_vio_space
    A1_ind
    intercept
    str_thol
    alpha
    params
    function_hatmatrix
    ncores
    function(colnames.cluster) {
      tryCatch_WEM(tsci_fit(
        df_treatment = df_treatment,
        Y = Y,
        D = D,
        Z = Z,
        X = X,
        list_vio_space = list_vio_space,
        intercept = intercept,
        str_thol = str_thol,
        split_prop = split_prop,
        alpha = alpha,
        params = params,
        function_hatmatrix = function_hatmatrix
      ), tsci_fit_NA_return(Q = list_vio_space$Q))}
  })

  if (do_parallel) {
    if (parallel == "multicore") {
      list_outputs <- parallel::mclapply(seq_len(nsplits), tsci_parallel, mc.cores = ncores)
    } else if (parallel == "snow") {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncores))
        parallel::clusterExport(cl, varlist = getNamespaceExports("TSML"))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        list_outputs <- parallel::parLapply(cl, seq_len(nsplits), tsci_parallel)
        parallel::stopCluster(cl)
        cl <- NULL # overwrite object which is responsible for the connection
      } else list_outputs <- parallel::parLapply(cl, seq_len(nsplits), tsci_parallel)
    }
  } else list_outputs <- lapply(seq_len(nsplits), tsci_parallel)

  check_list_outputs <- check_output(list_outputs = list_outputs, ind_start = 1)

  if (check_list_outputs$prop_na > 0.25) {
    stop(paste0("In more then 25% of the sample splits the output statistics could not be calculated.",
               check_list_outputs$error_string), call. = FALSE)
  }

  if (check_list_outputs$prop_na > 0) {
    nsplits_new <- ceiling(nsplits * 0.5)
    if (do_parallel) {
      if (parallel == "multicore") {
        list_outputs_new <- parallel::mclapply(seq_len(nsplits_new), tsci_parallel, mc.cores = ncores)
      } else if (parallel == "snow") {
        if (is.null(cl)) {
          cl <- parallel::makePSOCKcluster(rep("localhost", ncores))
          parallel::clusterExport(cl, varlist = getNamespaceExports("TSML"))
          if (RNGkind()[1L] == "L'Ecuyer-CMRG")
            parallel::clusterSetRNGStream(cl)
          list_outputs_new <- parallel::parLapply(cl, seq_len(nsplits_new), tsci_parallel)
          parallel::stopCluster(cl)
          cl <- NULL # overwrite object which is responsible for the connection
        } else list_outputs_new <- parallel::parLapply(cl, seq_len(nsplits_new), tsci_parallel)
      }
    } else list_outputs_new <- lapply(seq_len(nsplits_new), tsci_parallel)
    check_list_outputs_new <- check_output(list_outputs = list_outputs_new, ind_start = nsplits + 1)
    error_string <- paste0(check_list_outputs$error_string, check_list_outputs_new$error_string)
    if ((1 - check_list_outputs_new$prop_na) * nsplits_new < check_list_outputs$prop_na * nsplits) {
      stop(paste("Even after performing ",
                 nsplits + nsplits_new,
                 "sample splits there were still less than",
                 nsplits,
                 "for which the output statistics could be calculated.",
                 error_string), call. = FALSE)
    }
    pos_na <- which(check_list_outputs$ind_na)
    pos_stat <- which(!(check_list_outputs_new$ind_na))
    for(i in seq_len(length(pos_na))) {
      list_outputs[[pos_na[i]]] <- list_outputs_new[[pos_stat[i]]]
    }
    warning(paste0("In ",
                   sum(check_list_outputs$ind_na),
                   " of the ",
                   nrounds,
                   " data splits the output statistics could not be calculated.",
                   " Thus another ",
                   nrounds_new,
                   " data splits were performed.",
                   error_string), call. = FALSE)
  }

  aggregate_output(output_list = list_outputs, alpha = alpha, Q = list_vio_space$Q, mult_split_method = mult_split_method)
}
