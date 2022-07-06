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


  if (!is.null(object$mult_split_method)) {
    if (object$mult_split_method == "FWER") {
      sd_robust <- rep(".", length(object$Coef_robust))
      sd_all <- rep(".", length(object$Coef_all))
    } else {
      sd_robust <- object$sd_robust
      sd_all <- object$sd_all
    }
  } else {
    sd_robust <- object$sd_robust
    sd_all <- object$sd_all
  }

  coef_df <-
    data.frame(c(object$Coef_robust, object$Coef_all),
               c(sd_robust, sd_all),
               c(object$CI_robust[1, ], object$CI_all[1, ]),
               c(object$CI_robust[2, ], object$CI_all[2, ]),
               c(object$pval_robust, object$pval_all))
  colnames(coef_df) <- c("Estimate", "Std_Error", paste(100 * object$alpha/2, "%"), paste(100*(1 - object$alpha/2), "%"), "Pr(>|t|)")
  rownames(coef_df) <- names(c(object$Coef_robust, object$Coef_all))

  TreatmentModel_df <- data.frame("Estimation_Method" = object$FirstStage_model,
                                  "Residual_Standard_Error" = object$FirstStage_rse,
                                  "R_Squared" = object$FirstStage_Rsquared)

  OutcomeModel_df <- data.frame("Residual_Standard_Error" = object$SecondStage_rse,
                                "R_Squared" = object$SecondStage_Rsquared)
  row.names(OutcomeModel_df) <- names(object$SecondStage_rse)

  ViolationSpace_DF <- data.frame("q_comp" = object$q_comp,
                                  "q_robust" = object$q_robust,
                                  "Qmax" = object$Qmax)
  row.names(ViolationSpace_DF) <- names(object$Qmax)

  IVStrength_df <- data.frame("IV_Strength" = object$iv_str,
                              "IV_Threshold" = object$iv_thol)
  row.names(IVStrength_df) <- names(object$iv_str)


  output_list <- list(coefficients = coef_df,
                      iv_strength = IVStrength_df,
                      viospace_selection = ViolationSpace_DF,
                      treatment_model = TreatmentModel_df,
                      outcome_model = OutcomeModel_df,
                      sample_size_A1 = object$n_A1,
                      sample_size_A2 = object$n_A2,
                      n_splits = object$nsplits,
                      mult_split_method = object$mult_split_method,
                      alpha = object$alpha)
  class(output_list) <- "summary.tsci"
  return(output_list)
}
