#' Performs violation space selection and estimates the effect of the treatment on the outcome.
#'
#' @param Y the outcome vector.
#' @param D the treatment vector.
#' @param Cov_aug a matrix containing the violation space and the baseline covariates.
#' @param A1_ind indices of samples in A1.
#' @param weight n.A1 by n.A1 hat matrix of the first stage.
#' @param Q number of violation spaces to test, including no violation.
#' @param rm_ind the indices to remove to identify different violation space.
#' @param intercept logical, if \code{TRUE} an intercept is included in the outcome model.
#' @param layer logical, if \code{TRUE} violation space selection is performed.
#' @param str_thol numeric, the minimal value of the threshold of IV strength test.
#' @param alpha numeric, the significant level used for the CI intervals.
#' @param method the method that was used to get the hat matrix of the first stage.
#' @param alternative_bias_adjustment logical, if \code{TRUE} a slightly
#' different bias adjustment is used.
#'
#' @return a list containing the following the following components: \tabular{ll}{
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
#' }
#' @noRd
#'
#' @examples
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n)
#' D <- Z + 2 * Z^2 + X + rnorm(n)
#' Y <- D + Z + X + rnorm(n)
#' Q <- 4
#' vio_space <- matrix(NA,n,0)
#' for (q in 1:(Q-1)) {
#'   vio_space <- cbind(Z^q,vio_space)
#' }
#' rm_ind <- rep(list(NA),Q-1)
#' for (i in 1:(Q-1)) {
#'   rm_ind[[i]] = 1:(Q-i)
#' }
#' Data <- data.frame("D" = D, "Z" = Z, "X" = X)
#' Data_A1 <- Data[c(1:50), ]
#' Data_A2 <- Data[-c(1:50), ]
#' params_xgboost <- list("nrounds" = 10, "eta" = 0.3, "max_depth" = 2, "subsample" = 1,
#' "colsample_bytree" = 1)
#' hatmatrix_xgboost <- get_l2boost_hatmatrix(Data_A1 = Data_A1, Data_A2 = Data_A2,
#' params = params_xgboost, xgboost = TRUE)
#' Cov_aug <- cbind(vio_space, X)
#' tsci_secondstage_selection(Y = Y, D = D, Cov_aug = Cov_aug, A1_ind = c(1:50),
#'   weight = hatmatrix_xgboost$weight, Q = Q, rm_ind = rm_ind, intercept = TRUE,
#'   layer = TRUE, str_thol = 20, alpha = 0.05, method = "BO")
tsci_secondstage_selection_alternative <- function(Y,
                                       D,
                                       Cov_aug,
                                       A1_ind,
                                       weight,
                                       Q,
                                       rm_ind,
                                       intercept,
                                       layer,
                                       str_thol,
                                       alpha,
                                       method,
                                       alternative_bias_adjustment = FALSE) {
  Y <- as.matrix(Y)
  D <- as.matrix(D)
  Cov_aug <- as.matrix(Cov_aug)
  # constants
  n <- length(Y)
  n_A1 <- length(A1_ind)
  r_aug <- NCOL(Cov_aug)
  Y_A1 <- Y[A1_ind]
  D_A1 <- D[A1_ind]
  Cov_aug_A1 <- Cov_aug[A1_ind, ]
  # compute the representations
  Y_rep <- as.matrix(weight %*% Y_A1)
  D_rep <- as.matrix(weight %*% D_A1)
  Cov_rep <- as.matrix(weight %*% Cov_aug_A1)
  # the noise of treatment model
  delta_hat <- D_A1 - D_rep
  SigmaSqD <- mean(delta_hat^2)
  # save estimates for selection part
  if (method == "RF") {
    names <- c(paste("RF-q", seq(0, Q - 1), sep = ""),
               paste("RF-Cor-q", seq(0, Q - 1), sep = ""))
  } else if (method == "BO") {
    names <- c(paste("BO-q", seq(0, Q - 1), sep = ""),
               paste("BO-Cor-q", seq(0, Q - 1), sep = ""))
  } else if (method == "NN") {
    names <- c(paste("NN-q", seq(0, Q - 1), sep = ""),
               paste("NN-Cor-q", seq(0, Q - 1), sep = ""))
  } else {
    stop("No valid method selected. Choose either 'BO', 'NN' or 'RF'.")
  }

  Coef_vec <- sd_vec <- rep(NA, 2 * Q)
  names(Coef_vec) <- names(sd_vec) <- names
  # IV strength test and signal strength test
  iv_str <- iv_thol <- rep(NA, Q)
  names(iv_str) <- names(iv_thol) <- paste("q", seq(0, Q - 1), sep = "")
  # the noise of outcome model
  eps_hat <- rep(list(NA), Q)
  # the numerator of variance
  explained_iv <- rep(NA, Q)
  names(explained_iv) <- paste("q", seq(0, Q - 1), sep = "")
  trace_T <- explained_iv
  SigmaSqY <- trace_T
  # (I - P_V,W)
  P_perp <- rep(list(NA), Q)
  # t(I - Omega)%*%(I - P_V,W)
  A <- rep(list(NA), Q)

  ### fixed violation space, compute necessary inputs of selection part
  # save D_resid for the computation of H and z_alpha
  D_resid <- RSS_vec <- rep(list(NA), Q)
  for (index in seq_len(Q)) {
    if (index == Q) {
      if (intercept) {
        reg_ml <- stats::lm(Y_rep ~ D_rep + Cov_rep)
        betaHat <- stats::coef(reg_ml)[2]
      } else {
        reg_ml <- stats::lm(Y_rep ~ D_rep + Cov_rep - 1)
        betaHat <- stats::coef(reg_ml)[1]
      }
      Coef_vec[index] <- betaHat
      eps_hat[[index]] <- stats::resid(stats::lm(Y_A1 - D_A1 * betaHat ~ Cov_aug_A1))
      if (alternative_bias_adjustment == T) {
        P_perp[[index]] <- diag(n_A1) - Cov_aug_A1 %*%
          chol2inv(chol(t(Cov_aug_A1) %*% Cov_aug_A1)) %*% t(Cov_aug_A1)
        A[[index]] <- t(diag(n_A1) - weight) %*% P_perp[[index]]
      }
      SigmaSqY[index] <- mean(eps_hat[[index]]^2)
      stat_outputs <- tsci_secondstage_stats(D_rep,
                                             Cov_rep,
                                             weight,
                                             n,
                                             eps_hat[[index]],
                                             delta_hat,
                                             str_thol = str_thol)
    } else {
      if (intercept) {
        reg_ml <- stats::lm(Y_rep ~ D_rep + Cov_rep[, -rm_ind[[index]]])
        betaHat <- stats::coef(reg_ml)[2]
      } else {
        reg_ml <- stats::lm(Y_rep ~ D_rep + Cov_rep[, -rm_ind[[index]]] - 1)
        betaHat <- stats::coef(reg_ml)[1]
      }
      Coef_vec[index] <- betaHat
      eps_hat[[index]] <-
        stats::resid(stats::lm(Y_A1 - D_A1 * betaHat ~ Cov_aug_A1[, -rm_ind[[index]]]))
      if (alternative_bias_adjustment == T) {
        P_perp[[index]] <- diag(n_A1) - Cov_aug_A1 %*%
          chol2inv(chol(t(Cov_aug_A1) %*% Cov_aug_A1)) %*% t(Cov_aug_A1)
        A[[index]] <- t(diag(n_A1) - weight) %*% P_perp[[index]]
      }
      SigmaSqY[index] <- mean(eps_hat[[index]]^2)
      stat_outputs <- tsci_secondstage_stats(D_rep,
                                             Cov_rep[, -rm_ind[[index]]],
                                             weight,
                                             n,
                                             eps_hat[[index]],
                                             delta_hat,
                                             str_thol = str_thol)
    }
    # the statistics
    sd_vec[index] <- stat_outputs$sd
    sd_vec[index + Q] <- stat_outputs$sd
    iv_str[index] <- stat_outputs$iv_str
    iv_thol[index] <- stat_outputs$iv_thol
    explained_iv[index] <- stat_outputs$explained_iv
    D_resid[[index]] <- stat_outputs$D_resid
    trace_T[index] <- stat_outputs$trace_T
    RSS_vec[[index]] <- stat_outputs$RSS_V
  }
  # Residual sum of squares of D_rep~Cov_rep
  D_RSS <- iv_str * SigmaSqD

  # violation space selection
  # all of the q below are from 0 to (Q-1), so use q+1 to index the columns
  # Comparison and robust estimators
  Coef_robust <- sd_robust <- rep(NA, 4)
  if (method == "RF") {
    names(Coef_robust) <- names(sd_robust) <-
      c("RF-comp", "RF-Cor-comp", "RF-robust", "RF-Cor-robust")
  } else if (method == "BO") {
    names(Coef_robust) <- names(sd_robust) <-
      c("BO-comp", "BO-Cor-comp", "BO-robust", "BO-Cor-robust")
  } else if (method == "NN") {
    names(Coef_robust) <- names(sd_robust) <-
      c("NN-comp", "NN-Cor-comp", "NN-robust", "NN-Cor-robust")
  } else {
    stop("No valid method selected. Choose either 'BO', 'NN' or 'RF'.")
  }

  ivtest_vec <- (iv_str >= iv_thol)
  run_OLS <- weak_iv <- FALSE
  if (sum(ivtest_vec) == 0) {
    warning("Weak IV, even if the IV is assumed to be valid; run OLS") # stop, output results of OLS
    run_OLS <- TRUE
    Qmax <- 1
  } else {
    Qmax <- sum(ivtest_vec) - 1
    if (Qmax == 0) {
      warning("Weak IV, if the IV is invalid. We still test the IV invalidity.")
      Qmax <- 1
      weak_iv <- TRUE
    }
  }

  # calculate Covariance at Qmax and compute bias correction
  eps_Qmax <- eps_hat[[Qmax + 1]]
  Coef_Qmax <- rep(NA, Q)
  if (alternative_bias_adjustment == TRUE) weightSq <- t(weight) %*% weight
  for (i in seq_len(Q)) {
    if (alternative_bias_adjustment == TRUE) {
      adjustment_Qmax <- 1 / max(0.5, diag(A[[Q]]) -
                                   sum(RSS_vec[[i]] * diag(A[[Q]]) * delta_hat^2 /
                                         (1 - 2 * diag(weight) + diag(weightSq))) /
                                   D_RSS[i])
      adjustment_q <- 1 / max(0.5, diag(A[[i]]) -
                                sum(RSS_vec[[i]] * diag(A[[i]]) * delta_hat^2 /
                                      (1 - 2 * diag(weight) + diag(weightSq))) /
                                D_RSS[i])
      Coef_Qmax[i] <- Coef_vec[i] -
        sum(RSS_vec[[i]] * delta_hat * eps_Qmax / adjustment_Qmax) / D_RSS[i]
      Coef_vec[i + Q] <- Coef_vec[i] -
        sum(RSS_vec[[i]] * delta_hat * eps_hat[[i]] / adjustment_q) / D_RSS[i]
    } else {
      Coef_Qmax[i] <- Coef_vec[i] -
        sum(RSS_vec[[i]] * delta_hat * eps_Qmax) / D_RSS[i]
      Coef_vec[i + Q] <- Coef_vec[i] -
        sum(RSS_vec[[i]] * delta_hat * eps_hat[[i]]) / D_RSS[i]
    }
  }
  sd_vec[-seq_len(Q)] <- sd_vec[seq_len(Q)]

  ### Selection
  # define comparison matrix
  H <- beta_diff <- matrix(0, Qmax, Qmax)
  # compute H matrix
  for (q1 in seq(0, Qmax - 1)) {
    for (q2 in seq(q1 + 1, Qmax)) {
      H[q1 + 1, q2] <- as.numeric(sum((weight %*% D_resid[[q1 + 1]])^2 * eps_Qmax^2) /
                                    (D_RSS[q1 + 1]^2) +
                                    sum((weight %*% D_resid[[q2 + 1]])^2 * eps_Qmax^2) /
                                    (D_RSS[q2 + 1]^2) -
                                    2 * sum(eps_Qmax^2 * (weight %*% D_resid[[q1 + 1]]) *
                                              (weight %*% D_resid[[q2 + 1]])) / (D_RSS[q1 + 1] * D_RSS[q2 + 1]))
    }
  }
  # compute beta difference matrix, use Qmax
  for (q in seq(0, Qmax - 1)) {
    beta_diff[q + 1, seq(q + 1, Qmax)] <-
      abs(Coef_Qmax[q + 1] - Coef_Qmax[seq(q + 2, Qmax + 1)]) # use bias-corrected estimator
  }
  # bootstrap for the quantile of the differences
  max_val <- rep(NA, 300)
  eps_Qmax_cent <- eps_Qmax - mean(eps_Qmax)
  for (i in seq_len(300)) {
    diff_mat <- matrix(0, Qmax, Qmax)
    eps <- rep(NA, n_A1)
    for (j in seq_len(n_A1)) {
      U_j <- stats::rnorm(1)
      eps[j] <- eps_Qmax_cent[j] * U_j
    }
    eps_rep <- weight %*% eps
    for (q1 in seq(0, Qmax - 1)) {
      for (q2 in seq(q1 + 1, Qmax)) {
        diff_mat[q1 + 1, q2] <- sum(D_resid[[q2 + 1]] * eps_rep) /
          (D_RSS[q2 + 1]) - sum(D_resid[[q1 + 1]] * eps_rep) / (D_RSS[q1 + 1])
      }
    }
    diff_mat <- abs(diff_mat) / sqrt(H)
    max_val[i] <- max(diff_mat, na.rm = TRUE)
  }
  z_alpha <- 1.01 * stats::quantile(max_val, 0.975)
  diff_thol <- z_alpha * sqrt(H)
  # comparison matrix
  C_alpha <- ifelse(beta_diff <= diff_thol, 0, 1)
  # layer selection or not
  if (layer == TRUE) {
    # a vector indicating the selection of each layer
    sel_vec <- apply(C_alpha, 1, sum)
    if (all(sel_vec != 0)) {
      q_comp <- Qmax
    } else {
      q_comp <- min(which(sel_vec == 0)) - 1
    }
  } else {
    ### What if Q = 3(q2) and Qmax = 1?
    sel_val <- C_alpha[1, Qmax]
    if (sel_val == 1) {
      q_comp <- Qmax
    } else {
      q_comp <- 0
    }
  }

  ### invalidity of TSLS
  if (q_comp >= 1) {
    invalidity <- 1
  } else {
    invalidity <- 0
  }
  q_robust <- min(q_comp + 1, Qmax)
  Coef_robust[1] <- Coef_vec[q_comp + 1]
  Coef_robust[2] <- Coef_vec[q_comp + Q + 1]
  Coef_robust[3] <- Coef_vec[q_robust + 1]
  Coef_robust[4] <- Coef_vec[q_robust + Q + 1]
  sd_robust[1] <- sd_vec[q_comp + 1]
  sd_robust[2] <- sd_vec[q_comp + Q + 1]
  sd_robust[3] <- sd_vec[q_robust + 1]
  sd_robust[4] <- sd_vec[q_robust + Q + 1]
  CI_robust <- rbind(Coef_robust + stats::qnorm(alpha / 2) * sd_robust,
                     Coef_robust + stats::qnorm(1 - alpha / 2) * sd_robust)
  rownames(CI_robust) <- c("lower", "upper")
  returnList <- list(
    Coef_vec = Coef_vec,
    sd_vec = sd_vec,
    Coef_robust = Coef_robust,
    sd_robust = sd_robust,
    CI_robust = CI_robust,
    iv_str = iv_str, iv_thol = iv_thol,
    SigmaSqD = SigmaSqD,
    SigmaSqY = SigmaSqY,
    SigmaSqY_Qmax = mean(eps_Qmax^2),
    trace_T = trace_T,
    explained_iv = explained_iv,
    Qmax = Qmax,
    q_comp = q_comp, q_robust = q_robust,
    invalidity = invalidity,
    run_OLS = run_OLS,
    weak_iv = weak_iv
  )
  returnList
}
