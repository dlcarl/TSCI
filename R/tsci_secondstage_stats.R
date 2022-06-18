#' @title Statistics Calculation for violation space selection
#' @description Calculate the statistics needed for violation space selection
#'
#' @param D_rep transformed treatment of dimension n_A1 by 1 corresponding to a violation space
#' @param Cov_rep transformed augmented covariates (vio_space, X) of dimension n_A1 by (p+ncol(vio_space)) corresponding to a violation space
#' @param weight n_A1 by n_A1 weight matrix
#' @param eps_hat residuals in the outcome model
#' @param delta_hat residuals in Random Forest corresponding to samples in A1
#' @param str_thol the minimal value of the threshold of IV strength test
#'
#' @return:
#'     \item{\code{sd}}{estimated standard error for betaHat}
#'     \item{\code{D_resid}}{residuals for D_rep ~ Cov_rep}
#'     \item{\code{iv_str}}{estimated IV strength}
#'     \item{\code{iv_thol}}{threshold for IV strength test}
#'     \item{\code{diag_M}}{diagonal of matrix M_{RF}(V) corresponding to a violation space}
#' @noRd
#'
#' @importFrom stats resid lm rnorm quantile
#'
tsci_secondstage_stats <- function(D_rep, Cov_rep, weight, eps_hat, delta_hat, str_thol) {
  n_A1 <- length(D_rep)
  # compute the trace of M_{RF}(V)
  # the trace of M_{RF}(V) matrix can be computed as RSS of regressing each column of weight matrix on Cov.rep
  SigmaSqD <- mean(delta_hat^2)
  diag_M <- rep(NA, n_A1)
  for (j in seq_len(n_A1)) {
    diag_M[j] <- sum(resid(lm(weight[, j] ~ Cov_rep))^2)
  }
  trace_M <- sum(diag_M)
  D_rep2 <- weight %*% D_rep
  D_resid <- resid(lm(D_rep ~ Cov_rep))
  D_RSS <- sum(D_resid^2)
  iv_str <- D_RSS / SigmaSqD
  # this is the numerator of the variance of betaHat
  sd <- sqrt(sum(eps_hat^2 * (weight %*% D_resid)^2)) / D_RSS

  # bootstrap for the threshold of IV strength test
  boot_vec <- rep(NA, 300)
  delta_cent <- delta_hat - mean(delta_hat)
  for (i in seq_len(300)) {
    delta <- delta_cent * rnorm(n_A1)

    delta_rep <- weight %*% delta
    delta_resid <- resid(lm(as.matrix(delta_rep) ~ Cov_rep))
    boot_vec[i] <- sum(delta_resid^2) + 2 * sum(D_rep2 * delta_resid)
  }
  iv_thol <- quantile(boot_vec, 0.975) /
    SigmaSqD + max(2 * trace_M, str_thol)
  returnList <- list(
    sd = sd,
    D_resid = D_resid,
    iv_str = iv_str,
    iv_thol = iv_thol,
    diag_M = diag_M
  )
  returnList
}
