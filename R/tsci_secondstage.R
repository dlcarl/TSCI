#' @title Violation Space Selection
#' @description Selects violation space. It allows for user specific hat matrices (\code{weight}). Otherwise call \code{tsci_boosting} or \code{tsci_forest}.
#'
#' @param Y outcome with dimension n by 1
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param vio_space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation space selection will be performed according to provided violation space, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case
#' @param A1_ind the indices of samples in A1
#' @param weight n_A1 by n_A1 weight matrix
#' @param intercept logic, to include intercept in the outcome model or not
#' @param str_thol the minimal value of the threshold of IV strength test
#' @param alpha the significance level
#'
#' @return
#' \describe{
#'     \item{\code{Coef_all}}{a series of point estimators of treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{sd_all}}{standard errors of Coef_all}
#'     \item{\code{CI_all}}{confidence intervals for the treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{Coef_robust}}{the point estimators corresponding to the violation space selected by the robust comparison}
#'     \item{\code{sd_robust}}{the standard errors of Coef_robust}
#'     \item{\code{CI_robust}}{confidence intervals for the treatment effect with the violation space selected by the robust comparison}
#'     \item{\code{iv_str}}{IV strength corresponding to different violation spaces}
#'     \item{\code{iv_thol}}{the threshold of IV strength test corresponding to different violation spaces}
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the first violation space and run TSRF only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q_hat}}{the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' }
#' @export
#' @examples
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n)
#' D <- Z + X + rnorm(n)
#' Y <- D + X + rnorm(n) # no violation
#' A1_ind <- seq_len(n)
#' weight <- X %*% chol(chol2inv(t(X) %*% X)) %*% t(X)
#' tsci_secondstage(Y = Y, D = D, Z = Z, X = X, A1_ind = A1_ind, weight = weight)
#'
#' @importFrom stats coef lm qnorm quantile resid rnorm
tsci_secondstage <- function(Y,
                             D,
                             Z,
                             X,
                             vio_space = NULL,
                             A1_ind,
                             weight,
                             intercept = TRUE,
                             str_thol = 10,
                             alpha = 0.05) {
  Y <- as.matrix(Y)
  D <- as.matrix(D)
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  n <- length(Y)
  n_A1 <- length(A1_ind)

  if (nrow(weight) != ncol(weight)) {
    stop("Transformation matrix must be a square matrix")
  }
  if (nrow(weight) != n_A1) {
    stop("The samples to construct transformation matrix must be the same as samples in A1")
  }

  Y_A1 <- Y[A1_ind]
  D_A1 <- D[A1_ind]

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

  # define the augmentation of covariates as combination of violation space and baseline covariates
  Cov_aug <- cbind(vio_space, X)
  Cov_aug_A1 <- Cov_aug[A1_ind, ]

  # compute the representations
  Y_rep <- as.matrix(weight %*% Y_A1)
  D_rep <- as.matrix(weight %*% D_A1)
  Cov_rep <- as.matrix(weight %*% Cov_aug_A1)

  # the noise of treatment model
  delta_hat <- D_A1 - D_rep
  SigmaSqD <- mean(delta_hat^2)
  # save estimates for selection part
  Coef_names <- c(
    paste("TSCI-q", seq(0, Q - 1), sep = ""),
    paste("TSCI-Cor-q", seq(0, Q - 1), sep = "")
  )
  Coef_all <- sd_all <- rep(NA, 2 * Q)
  names(Coef_all) <- names(sd_all) <- Coef_names

  # IV strength test
  iv_str <- iv_thol <- rep(NA, Q)
  names(iv_str) <- names(iv_thol) <- paste("q", seq(0, Q - 1), sep = "")
  # the noise of outcome model
  eps_hat <- rep(list(NA), Q)
  ### fixed violation space, compute necessary inputs of selection part
  D_resid <- diag_M_list <- rep(list(NA), Q)
  for (index in seq_len(Q)) {
    if (index == Q) {
      if (intercept) {
        reg_ml <- lm(Y_rep ~ D_rep + Cov_rep)
        betaHat <- coef(reg_ml)[2]
      } else {
        reg_ml <- lm(Y_rep ~ D_rep + Cov_rep - 1)
        betaHat <- coef(reg_ml)[1]
      }
      Coef_all[index] <- betaHat
      eps_hat[[index]] <- resid(lm(Y_A1 - D_A1 * betaHat ~ Cov_aug_A1))
      stat_outputs <- tsci_secondstage_stats(D_rep,
        Cov_rep,
        weight,
        n,
        eps_hat[[index]],
        delta_hat,
        str_thol = str_thol
      )
    } else {
      if (intercept) {
        reg_ml <- lm(Y_rep ~ D_rep + Cov_rep[, -rm_ind[[index]]])
        betaHat <- coef(reg_ml)[2]
      } else {
        reg_ml <- lm(Y_rep ~ D_rep + Cov_rep[, -rm_ind[[index]]] - 1)
        betaHat <- coef(reg_ml)[1]
      }
      Coef_all[index] <- betaHat
      eps_hat[[index]] <-
        resid(lm(Y_A1 - D_A1 * betaHat ~ Cov_aug_A1[, -rm_ind[[index]]]))
      stat_outputs <- tsci_secondstage_stats(D_rep,
        Cov_rep[, -rm_ind[[index]]],
        weight,
        n,
        eps_hat[[index]],
        delta_hat,
        str_thol = str_thol
      )
    }
    # the necessary statistics
    sd_all[index] <- stat_outputs$sd
    sd_all[index + Q] <- stat_outputs$sd
    iv_str[index] <- stat_outputs$iv_str
    iv_thol[index] <- stat_outputs$iv_thol
    D_resid[[index]] <- stat_outputs$D_resid
    diag_M_list[[index]] <- stat_outputs$diag_M
  }
  # Residual sum of squares of D_rep~Cov_rep
  D_RSS <- iv_str * SigmaSqD

  # violation space selection
  # all of the q below are from 0 to (Q-1), so use q+1 to index the columns
  # Comparison and robust estimators
  Coef_robust <- sd_robust <- rep(NA, 2)
  names(Coef_robust) <- names(sd_robust) <- c("TSCI-comp", "TSCI-robust")
  ivtest_vec <- (iv_str >= iv_thol)
  if (sum(ivtest_vec) == 0) {
    warning("Weak IV, even if the IV is assumed to be valid; run OLS")
    Qmax <- -1
  } else {
    Qmax <- sum(ivtest_vec) - 1
    if (Qmax == 0) {
      warning("Weak IV, if the IV is invalid. We still test the IV invalidity.")
    }
  }

  # Compute bias-corrected estimators
  for (i in 1:Q) {
    Coef_all[i + Q] <- Coef_all[i] - sum(diag_M_list[[i]] * delta_hat * eps_hat[[i]]) / D_RSS[i]
  }
  sd_all[-(1:Q)] <- sd_all[1:Q]

  # If IV test fails at q0 or q1, we do not need to do selection
  if (Qmax >= 1) { # selection
    eps_Qmax <- eps_hat[[Qmax + 1]]
    Coef_Qmax <- rep(NA, Q)
    for (i in 1:Q) {
      Coef_Qmax[i] <- Coef_all[i] - sum(diag_M_list[[i]] * delta_hat * eps_Qmax) / D_RSS[i]
    }

    ### Selection
    # define comparison matrix
    H <- beta_diff <- matrix(0, Qmax, Qmax)
    # compute H matrix
    for (q1 in 0:(Qmax - 1)) {
      for (q2 in (q1 + 1):(Qmax)) {
        H[q1 + 1, q2] <- as.numeric(sum((weight %*% D_resid[[q1 + 1]])^2 * eps_Qmax^2) /
                                      (D_RSS[q1 + 1]^2) + sum((weight %*% D_resid[[q2 + 1]])^2 * eps_Qmax^2) /
                                      (D_RSS[q2 + 1]^2) -
          2 * sum(eps_Qmax^2 * (weight %*% D_resid[[q1 + 1]]) * (weight %*% D_resid[[q2 + 1]])) /
            (D_RSS[q1 + 1] * D_RSS[q2 + 1]))
      }
    }
    # compute beta difference matrix, use Qmax
    for (q in 0:(Qmax - 1)) {
      beta_diff[q + 1, (q + 1):(Qmax)] <- abs(Coef_Qmax[q + 1] - Coef_Qmax[(q + 2):(Qmax + 1)]) # use bias-corrected estimator
    }
    # bootstrap for the quantile of the differences
    max_val <- rep(NA, 300)
    eps_Qmax_cent <- eps_Qmax - mean(eps_Qmax)
    for (i in 1:300) {
      diff_mat <- matrix(0, Qmax, Qmax)
      eps <- rep(NA, n_A1)
      for (j in 1:n_A1) {
        U_j <- rnorm(1)
        eps[j] <- eps_Qmax_cent[j] * U_j
      }
      eps_rep <- weight %*% eps
      for (q1 in 0:(Qmax - 1)) {
        for (q2 in (q1 + 1):(Qmax)) {
          diff_mat[q1 + 1, q2] <- sum(D_resid[[q2 + 1]] * eps_rep) /
            (D_RSS[q2 + 1]) - sum(D_resid[[q1 + 1]] * eps_rep) / (D_RSS[q1 + 1])
        }
      }
      diff_mat <- abs(diff_mat) / sqrt(H)
      max_val[i] <- max(diff_mat, na.rm = TRUE)
    }
    z_alpha <- 1.01 * quantile(max_val, 0.975)
    diff_thol <- z_alpha * sqrt(H)
    # comparison matrix
    C_alpha <- ifelse(beta_diff <= diff_thol, 0, 1)


    # vector indicating the selection of each layer
    sel_vec <- apply(C_alpha, 1, sum)
    if (all(sel_vec != 0)) {
      q_comp <- Qmax
    } else {
      q_comp <- min(which(sel_vec == 0)) - 1
    }
  } # selection

  ### invalidity of TSLS
  if (Qmax >= 1) {
    if (q_comp >= 1) {
      invalidity <- 1
    } else {
      invalidity <- 0
    }
  } else {
    invalidity <- 0
  }

  # OLS estimator
  OLS <- summary(lm(Y ~ D + X))$coefficients
  Coef_OLS <- OLS[2, 1]
  sd_OLS <- OLS[2, 2]
  CI_OLS <- c(Coef_OLS + qnorm(alpha / 2) * sd_OLS, Coef_OLS + qnorm(1 - alpha / 2) * sd_OLS)
  names(CI_OLS) <- c("lower", "upper")
  # add OLS to Coef_all
  Coef_all <- c(Coef_all[(Q + 1):(2 * Q)], Coef_OLS)
  names(Coef_all) <- c(Coef_names[(Q + 1):(2 * Q)], "OLS")
  sd_all <- c(sd_all[(Q + 1):(2 * Q)], sd_OLS)
  names(sd_all) <- c(Coef_names[(Q + 1):(2 * Q)], "OLS")

  # confidence intervals for all violation spaces
  CI_all <- rbind(Coef_all + qnorm(alpha / 2) * sd_all, Coef_all + qnorm(1 - alpha / 2) * sd_all)
  rownames(CI_all) <- c("lower", "upper")


  ### estimated violation space and corresponding estimator
  if (Qmax >= 1) {
    q_robust <- min(q_comp + 1, Qmax)
    q_hat <- c(q_comp, q_robust)
    names(q_hat) <- c("comp", "robust")
    Coef_robust[1] <- Coef_all[q_comp + 1]
    Coef_robust[2] <- Coef_all[q_robust + 1]
    sd_robust[1] <- sd_all[q_comp + 1]
    sd_robust[2] <- sd_all[q_robust + 1]
  }
  if (Qmax == 0) {
    q_hat <- c(0)
    names(q_hat) <- "q0"
    Coef_robust <- Coef_all[1]
    sd_robust <- sd_all[1]
  }
  if (Qmax == -1) {
    q_hat <- c(-1)
    names(q_hat) <- "OLS"
    Coef_robust <- Coef_all[Q + 1]
    sd_robust <- sd_all[Q + 1]
  }

  CI_robust <- rbind(
    Coef_robust + qnorm(alpha / 2) * sd_robust,
    Coef_robust + qnorm(1 - alpha / 2) * sd_robust
  )
  rownames(CI_robust) <- c("lower", "upper")

  # change format of Qmax, q_hat and q_robust. Only temporary. Will change structure of code.
  max_Q <- comp_q <- robust_q <- rep(0, Q + 1)
  names(max_Q) <- names(comp_q) <- names(robust_q) <- paste("q", seq(-1, Q - 1), sep = "")
  max_Q[Qmax + 2] <- 1
  if(length(q_hat) == 2) {
    comp_q[q_hat[1] + 2] <- 1
    robust_q[q_hat[2]+ 2] <- 1
  } else {
    comp_q[q_hat + 2] <- 1
    robust_q[q_hat + 2] <- 1
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
    Qmax = max_Q,
    q_comp = comp_q,
    q_robust = robust_q,
    invalidity = invalidity
  )
  returnList
}
