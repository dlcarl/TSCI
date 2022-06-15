#' tsci_fit_NA_return
#'
#' @param Q number of violation spaces to test.
#'
#' @return an object consisting of NA-entries only but of the same dimension as
# the return object of the function beta_crossfit
#' @noRd
#'
tsci_fit_NA_return <- function(Q) {
  list(
    Coef_all = rep(NA, Q + 1),
    sd_all = rep(NA, Q + 1),
    CI_all = matrix(NA, nrow = 2, ncol = Q + 1),
    Coef_robust = rep(NA, 2),
    sd_robust = rep(NA, 2),
    CI_robust = rep(NA, 2),
    iv_str = rep(NA, Q),
    iv_thol = rep(NA, Q),
    Qmax = NA,
    q_comp = NA,
    invalidity = NA
  )
}
