#' Computes the necessary statistics for TSCI.
#'
#' @param D_rep a vector denoting the treatment representation.
#' @param Cov_rep a matrix denoting the augmented instruments-covariates representation.
#' @param weight a matrix, the hat matrix obtained in the first stage.
#' @param n numeric, the full sample size.
#' @param eps_hat a vector containing the estimated noise of the outcome model.
#' @param delta_hat a vector containing the estimated noise of the treatment model.
#' @param str_thol numeric, the minimal value of the threshold of IV strength test.
#'
#' @return a list containing the following the following components: \tabular{ll}{
#' \code{sd} \tab the estimated standard error of the treatment parameter estimate. \cr
#' \tab \cr
#' \code{D_resid} \tab BLANK. \cr
#' \tab \cr
#' \code{iv_str} \tab the estimated instrumental variable strength. \cr
#' \tab \cr
#' \code{iv_thol} \tab the estimated necessary instrumental variable strength. \cr
#' \tab \cr
#' \code{explained_iv} \tab the numerator of the estimated standard error. \cr
#' \tab \cr
#' \code{trace_T} \tab the trace of the M matrix. \cr
#' \tab \cr
#' \code{RSS_V} \tab the diagonal entries of the M matrix. \cr
#' }
#'
#' @export
#'
#' @examples
#' NULL
tsci_secondstage_stats <- function(D_rep,
                                    Cov_rep,
                                    weight,
                                    n,
                                    eps_hat,
                                    delta_hat,
                                    str_thol) {
  n_A1 <- length(D_rep)
  r_aug <- NCOL(Cov_rep)
  # compute the trace of T(V)
  # the trace of T matrix can be computed as RSS of each column of Omega on Cov_rep
  SigmaSqD <- mean(delta_hat^2)
  RSS_V <- rep(NA, n_A1)
  for (j in seq_len(n_A1)) {
    RSS_V[j] <- sum(stats::resid(stats::lm(weight[, j] ~ Cov_rep))^2)
  }
  trace_T <- sum(RSS_V)
  D_rep2 <- weight %*% D_rep
  D_resid <- stats::resid(stats::lm(D_rep ~ Cov_rep))
  D_RSS <- sum(D_resid^2)
  iv_str <- D_RSS / SigmaSqD
  # this is the numerator of the variance of betaHat
  explained_iv <- as.numeric((t(D_resid) %*% t(weight)) %*% weight %*% D_resid)
  sd <- sqrt(sum(eps_hat^2 * (weight %*% D_resid)^2)) / D_RSS

  # bootstrap for the threshold of IV strength test
  boot_vec <- rep(NA, 300)
  delta_cent <- delta_hat - mean(delta_hat)
  for (i in seq_len(300)) {
    delta <- rep(NA, n_A1)
    for (j in seq_len(n_A1)) {
      U_j <- stats::rnorm(1)
      delta[j] <- delta_cent[j] * U_j
    }

    delta_rep <- weight %*% delta
    delta_resid <- stats::resid(stats::lm(as.matrix(delta_rep) ~ Cov_rep))
    boot_vec[i] <- sum(delta_resid^2) + 2 * sum(D_rep2 * delta_resid)
  }
  iv_thol <- stats::quantile(boot_vec, 0.975) /
    SigmaSqD + max(2 * trace_T, str_thol)
  returnList <- list(
    sd = sd,
    D_resid = D_resid,
    iv_str = iv_str,
    iv_thol = iv_thol,
    explained_iv = explained_iv,
    trace_T = trace_T,
    RSS_V = RSS_V
  )
  returnList
}
