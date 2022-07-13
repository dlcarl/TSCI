#' @noRd
tsci_fit_NA_return <- function(Q) {
  # this function returns an object consisting of NA-entries only but of the same dimension as
  # the return object of the function tsci_selection
  Coef_all <- sd_all <- pval_all <- rep(NA, Q + 1)
  names(Coef_all) <- names(sd_all) <- names(pval_all) <- c("OLS", paste("TSCI-Cor-q", seq(0, Q - 1), sep = ""))
  CI_all <- matrix(NA, nrow = 2, ncol = Q + 1)
  rownames(CI_all) <- c("lower", "upper")
  colnames(CI_all) <- names(Coef_all)
  Coef_robust <- sd_robust <- pval_robust <- rep(NA, 2)
  names(Coef_robust) <- names(sd_robust) <- names(pval_robust) <- c("TSCI-comp", "TSCI-robust")
  CI_robust <- matrix(NA, nrow = 2, ncol = 2)
  rownames(CI_robust) <- c("lower", "upper")
  colnames(CI_robust) <- names(Coef_robust)
  iv_str <- iv_thol <- rep(NA, Q)
  names(iv_str) <- names(iv_thol) <- paste("q", seq(0, Q - 1), sep = "")
  Qmax <- q_comp <- q_robust <- rep(NA, Q + 1)
  names(Qmax) <- names(q_comp) <- names(q_robust) <- c("OLS", paste("q", seq(0, Q - 1), sep = ""))
  invalidity <- rep(NA, 2)
  names(invalidity) <- c("valid", "invalid")
  FirstStage_rse <- FirstStage_Rsquared <- NA
  SecondStage_rse <- SecondStage_Rsquared <- rep(NA, Q + 1)
  names(SecondStage_rse) <- names(SecondStage_Rsquared) <- c("OLS", paste("TSCI-Cor-q", seq(0, Q - 1), sep = ""))


  list(
    Coef_all = Coef_all,
    sd_all = sd_all,
    pval_all = pval_all,
    CI_all = CI_all,
    Coef_robust = Coef_robust,
    sd_robust = sd_robust,
    pval_robust = pval_robust,
    CI_robust = CI_robust,
    iv_str = iv_str,
    iv_thol = iv_thol,
    Qmax = Qmax,
    q_comp = q_comp,
    q_robust = q_robust,
    invalidity = invalidity,
    FirstStage_rse = FirstStage_rse,
    FirstStage_Rsquared = FirstStage_Rsquared,
    SecondStage_rse = SecondStage_rse,
    SecondStage_Rsquared = SecondStage_Rsquared
  )
}
