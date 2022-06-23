
#' Violation Space Selection
#'
#' @param Y_A1 outcome with dimension n by 1
#' @param D_A1 treatment with dimension n by 1
#' @param X_A1 baseline covariates with dimension n by p
#' @param vio_space the \code{matrix} of the largest violation space
#' @param rm_ind a \code{list} containing the indices to remove to obtain the violation spaces to test for (including the null space).
#' @param Q the number of violation spaces (including the null space).
#' @param weight n_A1 by n_A1 weight matrix.
#' @param intercept logic, to include intercept in the outcome model or not.
#' @param str_thol the minimal value of the threshold of IV strength test.
#' @param alpha alpha the significance level.
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
#' @noRd
#'
#' @importFrom stats coef lm qnorm quantile resid rnorm
#'
tsci_selection <- function(Y_A1,
                           D_A1,
                           X_A1,
                           vio_space,
                           rm_ind,
                           Q,
                           weight,
                           intercept,
                           str_thol,
                           alpha) {
  Cov_aug_A1 <- cbind(vio_space, X_A1)
  Y_rep <- as.matrix(weight %*% Y_A1)
  D_rep <- as.matrix(weight %*% D_A1)
  Cov_rep <- as.matrix(weight %*% Cov_aug_A1)
  n_A1 <- NROW(Y_rep)

  # initialize output list
  output <- tsci_fit_NA_return(Q = Q)

  # the non bias corrected beta estimates
  Coef_all <- rep(NA, Q)

  # the noise of treatment model
  delta_hat <- D_A1 - D_rep
  SigmaSqD <- mean(delta_hat^2)

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
        eps_hat[[index]],
        delta_hat,
        str_thol = str_thol
      )
    }
    # the necessary statistics
    output$sd_all[index] <- stat_outputs$sd
    output$iv_str[index] <- stat_outputs$iv_str
    output$iv_thol[index] <- stat_outputs$iv_thol
    D_resid[[index]] <- stat_outputs$D_resid
    diag_M_list[[index]] <- stat_outputs$diag_M
  }
  # Residual sum of squares of D_rep~Cov_rep
  D_RSS <- output$iv_str * SigmaSqD

  # violation space selection
  # all of the q below are from 0 to (Q-1), so use q+1 to index the columns
  # Comparison and robust estimators
  Coef_robust <- sd_robust <- rep(NA, 2)
  names(Coef_robust) <- names(sd_robust) <- c("TSCI-comp", "TSCI-robust")
  ivtest_vec <- (output$iv_str >= output$iv_thol)
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
  for (i in seq_len(Q)) {
    output$Coef_all[i] <- Coef_all[i] - sum(diag_M_list[[i]] * delta_hat * eps_hat[[i]]) / D_RSS[i]
  }

  # If IV test fails at q0 or q1, we do not need to do selection
  if (Qmax >= 1) { # selection
    eps_Qmax <- eps_hat[[Qmax + 1]]
    Coef_Qmax <- rep(NA, Q)
    for (i in seq_len(Q)) {
      Coef_Qmax[i] <- Coef_all[i] - sum(diag_M_list[[i]] * delta_hat * eps_Qmax) / D_RSS[i]
    }

    ### Selection
    # define comparison matrix
    H <- beta_diff <- matrix(0, Qmax, Qmax)
    # compute H matrix
    for (q1 in seq_len(Qmax) - 1) {
      for (q2 in (q1 + 1):(Qmax)) {
        H[q1 + 1, q2] <- as.numeric(sum((weight %*% D_resid[[q1 + 1]])^2 * eps_Qmax^2) /
          (D_RSS[q1 + 1]^2) + sum((weight %*% D_resid[[q2 + 1]])^2 * eps_Qmax^2) /
            (D_RSS[q2 + 1]^2) -
          2 * sum(eps_Qmax^2 * (weight %*% D_resid[[q1 + 1]]) * (weight %*% D_resid[[q2 + 1]])) /
            (D_RSS[q1 + 1] * D_RSS[q2 + 1]))
      }
    }
    # compute beta difference matrix, use Qmax
    for (q in seq_len(Qmax) - 1) {
      beta_diff[q + 1, (q + 1):(Qmax)] <- abs(Coef_Qmax[q + 1] - Coef_Qmax[(q + 2):(Qmax + 1)]) # use bias-corrected estimator
    }
    # bootstrap for the quantile of the differences
    B <- 300
    eps_Qmax_cent <- as.vector(eps_Qmax - mean(eps_Qmax))
    eps_rep_matrix <- weight %*% (eps_Qmax_cent * matrix(rnorm(n_A1 * B), ncol = B))
    diff_mat <- matrix(0, Qmax, Qmax)
    max_val <-
      apply(eps_rep_matrix, 2,
            FUN = function(eps_rep) {
              for (q1 in seq_len(Qmax) - 1) {
                for (q2 in (q1 + 1):(Qmax)) {
                  diff_mat[q1 + 1, q2] <- sum(D_resid[[q2 + 1]] * eps_rep) /
                    (D_RSS[q2 + 1]) - sum(D_resid[[q1 + 1]] * eps_rep) / (D_RSS[q1 + 1])
                }
              }
              diff_mat <- abs(diff_mat) / sqrt(H)
              max(diff_mat, na.rm = TRUE)
            })
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
      output$invalidity[1] <- 0
      output$invalidity[2] <- 1
    } else {
      output$invalidity[1] <- 1
      output$invalidity[2] <- 0
    }
  } else {
    output$invalidity[1] <- 1
    output$invalidity[2] <- 0
  }

  # OLS estimator
  OLS <- summary(lm(Y ~ D + X))$coefficients
  Coef_OLS <- OLS[2, 1]
  sd_OLS <- OLS[2, 2]
  # add OLS to Coef_all
  output$Coef_all[Q + 1] <- Coef_OLS
  output$sd_all[Q + 1] <- sd_OLS

  # confidence intervals and p values for all violation spaces
  output$CI_all[] <-
    rbind(output$Coef_all + qnorm(alpha / 2) * output$sd_all,
          output$Coef_all + qnorm(1 - alpha / 2) * output$sd_all)
  output$pval_all[] <-
    sapply(seq_len(length(output$Coef_all)),
      FUN = function(i) p_val(Coef = output$CI_all[i], SE = output$sd_all[i], beta_test = 0)
    )

  ### estimated violation space and corresponding estimator
  output$Qmax[] <- 0
  output$q_comp[] <- 0
  output$q_robust[] <- 0
  if (Qmax >= 1) {
    q_robust <- min(q_comp + 1, Qmax)
    output$Coef_robust[1] <- output$Coef_all[q_comp + 1]
    output$Coef_robust[2] <- output$Coef_all[q_robust + 1]
    output$sd_robust[1] <- output$sd_all[q_comp + 1]
    output$sd_robust[2] <- output$sd_all[q_robust + 1]
  } else if (Qmax == 0) {
    q_comp <- 0
    q_robust <- 0
    output$Coef_robust[1] <- output$Coef_all[q_comp + 1]
    output$Coef_robust[2] <- output$Coef_all[q_robust + 1]
    output$sd_robust[1] <- output$sd_all[q_comp + 1]
    output$sd_robust[2] <- output$sd_all[q_robust + 1]
  } else if (Qmax == -1) {
    q_comp <- -1
    q_robust <- -1
    output$Coef_robust[1] <- output$Coef_all[Q + 1]
    output$Coef_robust[2] <- output$Coef_all[Q + 1]
    output$sd_robust[1] <- output$sd_all[Q + 1]
    output$sd_robust[2] <- output$sd_all[Q + 1]
  }
  output$Qmax[Qmax + 2] <- 1
  output$q_comp[q_comp + 2] <- 1
  output$q_robust[q_robust + 2] <- 1

  output$CI_robust[] <- rbind(
    output$Coef_robust + qnorm(alpha / 2) * output$sd_robust,
    output$Coef_robust + qnorm(1 - alpha / 2) * output$sd_robust
  )
  output$pval_robust[] <-
    sapply(seq_len(length(output$Coef_robust)),
      FUN = function(i) p_val(Coef = output$CI_robust[i], SE = output$sd_robust[i], beta_test = 0)
    )

  output
}
