#' xxx
#'
#' @param object xxxx
#' @param ... xxx
#'
#' @return xxx
#' @export
summary.tsci <- function(object,
                         ...) {
  stopifnot(inherits(object, "tsci"))

  coef_robust <- object$Coef_robust
  sd_robust <- object$sd_robust
  pval_robust <- object$pval_robust
  ci_robust <- object$CI_robust

  coef_all <- object$Coef_all
  sd_all <- object$sd_all
  pval_all <- object$pval_all
  ci_all <- object$CI_all

  iv_str <- object$iv_str
  iv_tol <- object$iv_thol

  Qmax <- object$Qmax
  q_comp <- object$q_comp
  q_robust <- object$q_robust

  FirstStage_model <- object$FirstStage_model
  FirstStage_params <- object$FirstStage_params
  FirstStage_rse <- object$FirstStage_rse
  FirstStage_Rsquared <- object$FirstStage_Rsquared

  SecondStage_rse <- object$SecondStage_rse
  SecondStage_Rsquared <- object$SecondStage_Rsquared

  split_prop <- object$split_prop
  nsplits <- object$nsplits
  mult_split_method <- object$mult_split_method

  output_list <- list(coef_robust = coef_robust,
                      sd_robust = sd_robust,
                      pval_robust = pval_robust,
                      ci_robust = ci_robust,
                      coef_all = coef_all,
                      sd_all = sd_all,
                      pval_all = pval_all,
                      ci_all,
                      iv_str = iv_str,
                      iv_tol = iv_tol,
                      Qmax = Qmax,
                      q_comp = q_comp,
                      q_robust = q_robust,
                      FirstStage_model = FirstStage_model,
                      FirstStage_params = FirstStage_params,
                      FirstStage_rse = FirstStage_rse,
                      FirstStage_Rsquared = FirstStage_Rsquared,
                      SecondStage_rse = SecondStage_rse,
                      SecondStage_Rsquared = SecondStage_Rsquared,
                      split_prop = split_prop,
                      nsplits = nsplits,
                      mult_split_method = mult_split_method)
  class(output_list) <- "summary.tsci"
  return(output_list)
}
