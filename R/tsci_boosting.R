#' Calculates the hat matrix of the first stage.
#'
#' @param Y the outcome vector.
#' @param D the treatment vector.
#' @param Z the instrumental variable vector.
#' @param X the baseline covariates.
#' @param intercept logical, if \code{TRUE} an intercept is included in the outcome model.
#' @param vio_space list or matrix containing the violation space.
#' @param layer logical, if \code{TRUE} violation space selection is performed.
#' @param split_prop numeric, proportion of observations used to fit the outcome model.
#' @param nrounds numeric, hyperparameter of the boosting algorithm. Specifies the number of boosting iterations.
#' @param eta numeric, hyperparameter of the boosting algorithm. Specifies the learning rate.
#' @param max_depth numeric, hyperparameter of the boosting algorithm. Specifies the maximal depth of each tree.
#' @param subsample numeric, hyperparameter of the boosting algorithm. Specifies proportion of observations used to fit each tree.
#' @param colsample_bytree numeric, hyperparameter of the boosting algorithm. Specifies proportion of variables used to fit each tree.
#' @param early_stopping logical, hyperparameter of the boosting algorithm. If \code{TRUE} early stopping will be applied.
#' @param nfolds numeric, the number of folds for the k-fold cross validation.
#' @param l2boost_save logical, specifies if the fitted boosting model should be returned. Currently not implemented.
#' @param str_thol numeric, the minimal value of the threshold of IV strength test.
#' @param alpha numeric, specifies the desired significance level.
#' @param xgboost logical, specifies whether xgboost should be used to calculate the hat matrix.
#' @param alternative_bias_adjustment logical, if \code{TRUE} a slightly
#' different bias adjustment is used.
#' @param multi_splitting logical, if \code{TRUE} multi-splitting will be performed.
#' @param nsplits numeric, number of times the data will be split.
#' @param mult_split_method method to for inference if multi-splitting is performed. Either 'Chernozhukov' or 'Meinshausen'.
#'
#' @return a list containing at least the following the following components: \tabular{ll}{
#' \code{Coef_vec} \tab a named vector containing the coefficients for all tested violation spaces. \cr
#' \tab \cr
#' \code{sd_vec} \tab a named vector containing the standard errors of the coefficients
#' for all tested violation spaces. \cr
#' \tab \cr
#' \code{Coef_robust} \tab a named vector containing the coefficients of the selected violation space. \cr
#' \tab \cr
#' \code{sd_robust} \tab a named vector containing the standard errors of the coefficients
#' of the selected violation space. \cr
#' \tab \cr
#' \code{CI_robust} \tab a named vector containing the confidence intervals of the coefficients
#' of the selected violation space. \cr
#' \tab \cr
#' \code{iv_str} \tab a named vector containing the estimated instrumental variable strengths
#' of all tested violation spaces. \cr
#' \tab \cr
#' \code{iv_str} \tab a named vector containing the estimated necessary instrumental variable strengths
#' of all tested violation spaces. \cr
#' \code{SigmaSqD} \tab the estimated variance of the treatment model. \cr
#' \tab \cr
#' \code{SigmaSqY} \tab a named vector containing the estimated variances of the treatment model
#' of all tested violation spaces. \cr
#' \tab \cr
#' \code{SigmaSqY.Qmax} \tab the estimated variance of the treatment model
#' of the largest considered violation space. \cr
#' \tab \cr
#' \code{trace_T} \tab the trace of the M matrix. \cr
#' \tab \cr
#' \code{explained_iv} \tab a named vector containing the numerators of the estimated standard errors of the coefficients
#' of the selected violation space. \cr
#' \tab \cr
#' \code{Qmax} \tab the number of the largest considered violation space \cr
#' \tab \cr
#' \code{q_comp} \tab the number of the selected violation space \cr
#' \tab \cr
#' \code{q_robust} \tab the number of the selected violation space by the robust approach. \cr
#' \tab \cr
#' \code{invalidity} \tab logical, if \code{TRUE} the hypothesis of a valid instrument is rejected. \cr
#' \tab \cr
#' \code{run_OLS} \tab logical, if \code{TRUE} the instrument is weak even if it's assumed to be valid. \cr
#' \tab \cr
#' \code{run_OLS} \tab logical, if \code{TRUE} the instrument is weak for all tested non empty violation spaces. \cr
#' \tab \cr
#' \code{mse_cv} \tab the cross-validation mean squared error of the first stage.
#' }
#' @export
#'
#' @examples
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n)
#' D <- Z + 2 * Z^2 + X + rnorm(n)
#' Y <- D + Z + X + rnorm(n)
#' tsci_boosting(Y = Y, D = D, Z = Z, X = X)
tsci_boosting <- function(Y,
                          D,
                          Z,
                          X,
                          intercept = TRUE,
                          vio_space = NULL,
                          layer = TRUE,
                          split_prop = 2 / 3,
                          nrounds = NULL,
                          eta = NULL,
                          max_depth = NULL,
                          subsample = NULL,
                          colsample_bytree = NULL,
                          early_stopping = NULL,
                          nfolds = 10,
                          l2boost_save = TRUE,
                          str_thol = 20,
                          alpha = 0.05,
                          xgboost = TRUE,
                          alternative_bias_adjustment = FALSE,
                          multi_splitting = FALSE,
                          nsplits = NULL,
                          mult_split_method = "Meinshausen") {
  # reformat data
  Y <- as.matrix(Y)
  D <- as.matrix(D)
  Z <- as.matrix(Z)
  X <- as.matrix(X)

  # constants
  n <- NROW(X)
  p <- NCOL(X)

  # define defaults for the hyperparameters
  if (is.null(nrounds)) nrounds <- 50
  if (is.null(eta)) eta <- 0.3
  if (is.null(max_depth)) max_depth <- 6
  if (is.null(subsample)) subsample <- c(1)
  if (is.null(colsample_bytree)) colsample_bytree <- c(1)
  if (is.null(early_stopping)) early_stopping <- TRUE
  if (multi_splitting == TRUE) {
    split_prop <- 0.5
    if (is.null(nsplits)) nsplits <- 50
    if (!(mult_split_method %in% c("Meinshausen", "Chernozhukov"))) {
      stop("No valid multi-splitting inference method selected.
           Choose either 'Chernozhukov' or 'Meinshausen'.")
    }
  }

  # define the vio_space as polynomials if not specified
  if (is.null(vio_space)) {
    Q <- 4
    vio_space <- matrix(NA, n, 0)
    for (q in 1:(Q - 1)) {
      vio_space <- cbind(Z^q, vio_space)
    }
    # the indices to remove to identify violation space
    rm_ind <- rep(list(NA), Q - 1)
    for (i in 1:(Q - 1)) {
      rm_ind[[i]] <- 1:(Q - i)
    }
  }

  if (!is.null(vio_space)) {
    if (class(vio_space)[1] == "list") {
      vio_space <- lapply(vio_space, as.matrix)
      Q <- length(vio_space) + 1
      v_len <- sapply(vio_space, dim)[2, ]
      # the indices to remove to identify violation space
      rm_ind <- rep(list(NA), Q - 1)
      for (i in 1:(Q - 1)) {
        rm_ind[[i]] <- 1:sum(v_len[1:(Q - i)])
      }
      # merge the list of violation space to a matrix
      vio_space <- Reduce(cbind, vio_space)
    } else if (class(vio_space)[1] == "matrix") {
      Q <- ncol(vio_space) + 1
      rm_ind <- rep(list(NA), Q - 1)
      for (i in 1:(Q - 1)) {
        rm_ind[[i]] <- 1:(Q - i)
      }
    }
  }

  # define the augmentation of covariates,
  # which is the combination of violation space and baseline covariates
  Cov_aug <- cbind(vio_space, X)

  # Treatment model fitting
  W <- as.matrix(cbind(Z, X))
  D <- as.matrix(D)
  n <- NROW(W)
  p <- NCOL(W)
  Data <- data.frame(cbind(D, W))
  names(Data) <- c("D", paste("W", seq_len(p), sep = ""))
  # grid search
  params_grid <- expand.grid(
    nrounds = nrounds,
    eta = eta,
    max_depth = max_depth,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    early_stopping = early_stopping
  )

  # split the data into two parts A1 and A2
  # use A2 to train and use A1 to predict
  n_A1 <- round(split_prop * n)
  n_A2 <- n - n_A1
  A1_ind <- seq_len(n_A1)
  Data_A1 <- Data[A1_ind, ]
  Data_A2 <- Data[-A1_ind, ]

  # perform cross validation to select boosting hyperparameters
  treeboost_CV <- l2boost_CV(
    Data_A2 = Data_A2,
    params_grid = params_grid,
    nfolds = nfolds,
    xgboost = xgboost
  )


  if (multi_splitting == TRUE) {
    # initialize matrices
    Coef_matrix <- sd_matrix <- matrix(NA, nrow = nsplits, ncol = 2 * Q)
    iv_str_matrix <- iv_thol_matrix <- SigmaSqY_matrix <- trace_T_matrix <-
      explained_iv_matrix <- matrix(NA, nrow = nsplits, ncol = Q)
    Coef_robust_matrix <- sd_robust_matrix <- matrix(NA, nrow = nsplits, ncol = 4)
    SigmaSqD_matrix <- SigmaSqY_Qmax_matrix <- Qmax_matrix <- q_comp_matrix <-
      q_robust_matrix <- invalidity_matrix <- run_OLS_matrix <- weak_iv_matrix <-
      matrix(NA, nrow = nsplits, ncol = 1)

    for (s in seq_len(nsplits)) {
      # split the data into two parts A1 and A2
      # use A2 to train and use A1 to predict
      n_A1 <- round(split_prop * n)
      n_A2 <- n - n_A1
      A1_ind <- sample(seq_len(n), n_A1)
      Data_A1 <- Data[A1_ind, ]
      Data_A2 <- Data[-A1_ind, ]

      # refit model on whole training set A2 and calculate weight matrix for A1
      treeboost <- get_l2boost_hat_matrix(
        Data_A1 = Data_A1,
        Data_A2 = Data_A2,
        params = treeboost_CV$params,
        xgboost = xgboost
      )

      # Selection of violation space and estimation of treatment effect
      outputs <- tsci_second_stage_selection(Y,
        D,
        Cov_aug,
        A1_ind,
        weight = treeboost$weight,
        Q = Q,
        rm_ind = rm_ind,
        intercept = intercept,
        layer = layer,
        str_thol = str_thol,
        alpha = alpha,
        method = "BO",
        alternative_bias_adjustment = alternative_bias_adjustment
      )

      Coef_matrix[s, ] <- outputs$Coef_vec
      sd_matrix[s, ] <- outputs$sd_vec
      Coef_robust_matrix[s, ] <- outputs$Coef_robust
      sd_robust_matrix[s, ] <- outputs$sd_robust
      iv_str_matrix[s, ] <- outputs$iv_str
      iv_thol_matrix[s, ] <- outputs$iv_thol
      SigmaSqD_matrix[s, ] <- outputs$SigmaSqD
      SigmaSqY_matrix[s, ] <- outputs$SigmaSqY
      SigmaSqY_Qmax_matrix[s, ] <- outputs$SigmaSqY_Qmax
      trace_T_matrix[s, ] <- outputs$trace_T
      explained_iv_matrix[s, ] <- outputs$explained_iv
      Qmax_matrix[s, ] <- outputs$Qmax
      q_comp_matrix[s, ] <- outputs$q_comp
      q_robust_matrix[s, ] <- outputs$q_robust
      invalidity_matrix[s, ] <- outputs$invalidity
      run_OLS_matrix[s, ] <- outputs$run_OLS
      weak_iv_matrix[s, ] <- outputs$weak_iv
    }
    if (mult_split_method == "Meinshausen") {
      outputs$sd_vec[] <- NA
      outputs$sd_robust[] <- NA
      CI_Coef_vec <- sapply(seq_len(NCOL(Coef_matrix)),
        FUN = function(j) {
          beta <- Coef_matrix[, j]
          se <- sd_matrix[, j]
          # searches for the values of beta such that the median p.values are exactly half of alpha.
          # this may fail if the median of the estimated beta is not in the confidence interval
          lower <- stats::median(beta + -1 * stats::qnorm(1 - alpha / 4) * se)
          upper <- stats::median(beta + 1 * stats::qnorm(1 - alpha / 4) * se)
          return(c(lower, upper))
        }
      )
      rownames(CI_Coef_vec) <- c("lower", "upper")
      colnames(CI_Coef_vec) <- names(outputs$Coef_vec)
      outputs$CI_Coef_vec <- CI_Coef_vec
      CI_robust <- sapply(seq_len(NCOL(Coef_robust_matrix)),
        FUN = function(j) {
          beta <- Coef_robust_matrix[, j]
          se <- sd_robust_matrix[, j]
          # searches for the values of beta such that the median p.values are exactly half of alpha.
          # this may fail if the median of the estimated beta is not in the confidence interval
          lower <- stats::median(beta + -1 * stats::qnorm(1 - alpha / 4) * se)
          upper <- stats::median(beta + 1 * stats::qnorm(1 - alpha / 4) * se)
          return(c(lower, upper))
        }
      )
      outputs$CI_robust[] <- CI_robust
    } else if (mult_split_method == "Chernozhukov") {
      sd_vec <- sapply(seq_len(NCOL(Coef_matrix)),
        FUN = function(j) {
          stats::median(sqrt(sd_matrix[, j]^2 +
            (Coef_matrix[, j] - stats::median(Coef_matrix[, j]))^2))
        }
      )
      outputs$sd_vec[] <- sd_vec
      sd_robust <- sapply(seq_len(NCOL(Coef_robust_matrix)),
        FUN = function(j) {
          stats::median(sqrt(sd_robust_matrix[, j]^2 +
            (Coef_robust_matrix[, j] - stats::median(Coef_robust_matrix[, j]))^2))
        }
      )
      outputs$sd_robust[] <- sd_robust
      CI_robust <- sapply(seq_len(NCOL(Coef_robust_matrix)),
        FUN = function(j) {
          lower <- stats::median(Coef_robust_matrix[, j]) -
            stats::qnorm(1 - alpha / 2) * sd_robust[j]
          upper <- stats::median(Coef_robust_matrix[, j]) +
            stats::qnorm(1 - alpha / 2) * sd_robust[j]
          return(c(lower, upper))
        }
      )
      outputs$CI_robust[] <- CI_robust
    } else {
      outputs$sd_vec[] <- NA
      outputs$sd_robust[] <- NA
      outputs$CI_robust[] <- NA
    }
    outputs$Coef_vec[] <- apply(Coef_matrix, 2, stats::median)
    outputs$Coef_robust[] <- apply(Coef_robust_matrix, 2, stats::median)
    outputs$iv_str[] <- apply(iv_str_matrix, 2, stats::median)
    outputs$iv_thol[] <- apply(iv_thol_matrix, 2, stats::median)
    outputs$SigmaSqD[] <- apply(SigmaSqD_matrix, 2, stats::median)
    outputs$SigmaSqY[] <- apply(SigmaSqY_matrix, 2, stats::median)
    outputs$SigmaSqY_Qmax[] <- apply(SigmaSqY_Qmax_matrix, 2, stats::median)
    outputs$trace_T[] <- apply(trace_T_matrix, 2, stats::median)
    outputs$explained_iv[] <- apply(explained_iv_matrix, 2, stats::median)
    outputs$Qmax[] <- apply(Qmax_matrix, 2, stats::median)
    outputs$q_comp[] <- apply(q_comp_matrix, 2, stats::median)
    outputs$q_robust[] <- apply(q_robust_matrix, 2, stats::median)
    outputs$run_OLS[] <- apply(run_OLS_matrix, 2, stats::median)
    outputs$weak_iv[] <- apply(weak_iv_matrix, 2, stats::median)
  } else {
    # refit model on whole training set A2 and calculate weight matrix for A1
    treeboost <- get_l2boost_hat_matrix(
      Data_A1 = Data_A1,
      Data_A2 = Data_A2,
      params = treeboost_CV$params,
      xgboost = xgboost
    )

    # Selection of violation space and estimation of treatment effect
    outputs <- tsci_second_stage_selection(Y,
      D,
      Cov_aug,
      A1_ind,
      weight = treeboost$weight,
      Q = Q,
      rm_ind = rm_ind,
      intercept = intercept,
      layer = layer,
      str_thol = str_thol,
      alpha = alpha,
      method = "BO",
      alternative_bias_adjustment = alternative_bias_adjustment
    )
  }


  # Return output
  if (l2boost_save == TRUE) {
    outputs <- append(
      outputs,
      list(
        "mse_cv" = treeboost_CV$MSE_CV_A2,
        "l2boost" = treeboost$treeboost_A2
      )
    )
  } else {
    outputs <- append(
      outputs,
      list("mse_cv" = treeboost_CV$MSE_CV_A2)
    )
  }

  return(outputs)
}
