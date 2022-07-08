#' Title xxx
#'
#' @param df xxx
#' @param params_list xxx
#' @param gcv xxx
#' @param nfolds xxx
#'
#' @return xxx
#' @noRd
backfitting <- function(df,
                        params_list,
                        gcv,
                        nfolds) {
  max_iter <- 100
  tol <- 10^-6
  Y <- df[, 1]
  X <- df[, -1]
  n <- NROW(Y)
  p <- NCOL(X)

  mu <- mean(Y)
  g_X <- matrix(0, nrow = n, ncol = p)
  order_opt <- rep(NA, p)
  for (iter in seq_len(max_iter)) {
    g_old <- g_X
    for (j in seq_len(p)) {
      U_j <- Y - mu - rowSums(g_X[, -j])
      MSE_CV <- Inf
      order_j <- params_list[[j]]
      for (i in seq_len(length(order_j))) {
        A_j <- poly(X[, j], degree = order_j[i])
        if (gcv) {
          fit <- lm(U_j ~ A_j - 1)
          mse_cv <- mean((U_j - fit$fitted.values)^2) / (1 - NCOL(A_j) / NROW(A_j))^2
        } else {
          mse_cv <- cross_validate(Y = U_j, X = A_j, nfolds = nfolds)
        }
        if (mse_cv < MSE_CV) {
          order_opt[j] <- order_j[i]
          MSE_CV <- mse_cv
          g_X[, j] <- lm(U_j ~ A_j - 1)$fitted.values
        }
      }
    }
    change <- sqrt(sum((g_X - g_old)^2)) / sqrt(sum(g_old^2))
    if(change < tol) return(order_opt)
  }
  warning("backfitting did not converge.")
  return(order_opt)
}
