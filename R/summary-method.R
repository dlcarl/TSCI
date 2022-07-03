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
      sd_robust <- round(object$sd_robust, 5)
      sd_all <- round(object$sd_all, 5)
    }
  } else {
    sd_robust <- round(object$sd_robust, 5)
    sd_all <- round(object$sd_all, 5)
  }

  coef_df <-
    data.frame(round(c(object$Coef_robust, object$Coef_all), 5),
               c(sd_robust, sd_all),
               c(round(object$CI_robust[1, ], 5), round(object$CI_all[1, ], 5)),
               c(round(object$CI_robust[2, ], 5), round(object$CI_all[2, ], 5)),
               c(round(object$pval_robust, 3), round(object$pval_all, 3)))
  colnames(coef_df) <- c("Estimate", "Std. Error", paste(100 * object$alpha/2, "%"), paste(100*(1 - object$alpha/2), "%"), "Pr(>|t|)")
  rownames(coef_df) <- names(c(object$Coef_robust, object$Coef_all))

  TreatmentModel_df <- data.frame("Estimation Method" = object$FirstStage_model,
                                  "Residual Standard Error" = round(object$FirstStage_rse, 4),
                                  "R Squared" = round(object$FirstStage_Rsquared, 4))

  OutcomeModel_df <- data.frame("Residual Standard Error" = round(object$SecondStage_rse, 4),
                                "R Squared" = round(object$SecondStage_Rsquared, 4))
  row.names(OutcomeModel_df) <- names(object$SecondStage_rse)

  ViolationSpace_DF <- data.frame("q comp" = object$q_comp,
                                  "q robust" = object$q_robust,
                                  "Qmax" = object$Qmax)
  row.names(ViolationSpace_DF) <- names(object$Qmax)

  IVStrength_df <- data.frame("IV Strength" = object$iv_str,
                              "IV Tolerance" = object$iv_thol)
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
