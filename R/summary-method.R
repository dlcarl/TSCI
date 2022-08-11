#' Summarizing Two Stage Curvature Identification Fits
#'
#' @param object a object of class 'tsci'.
#' @param ... arguments to be passed to or from other methods.
#'
#' @return a object of class 'summary.tsci' containing the following elements:
#' \describe{
#'     \item{\code{coefficients}}{a data frame with columns for the estimated treatment coefficients, their standard errors, confidence intervals and (two-sided) p-values.}
#'     \item{\code{iv_strength}}{a data frame with columns containing the estimated instrumental variable strength and the estimated instrumental variable strength threshold.}
#'     \item{\code{viospace_selection}}{a data frame with columns for the frequency of violation space candidates selected by comparison,
#'     the robust method and as the largest violation space candidate for which the instrumental variable was considered to be strong enough.}
#'     \item{\code{treatment_model}}{a data frame with columns for the estimation method, the residual standard error and the R-squared of the fitted treatment model.}
#'     \item{\code{outcome_model}}{a data frame with columns for the residual standard error and the R-squared of the fitted outcome models.
#'     Except for OLS the values for R-squared are calculated keeping the treatment coefficient fixed. Thus, negative values are possible.
#'     If no intercept was included in the outcome model, values larger than 1 are also possible.}
#'     \item{\code{sample_size_A1}}{the number of observations in the subset used to fit the outcome model.}
#'     \item{\code{sample_size_A2}}{the number of observations in the subset used to train the parameters for fitting the treatment model.}
#'     \item{\code{n_splits}}{the number of sample splits performed.}
#'     \item{\code{mult_split_method}}{the method used to calculate the standard errors and p-values if \code{n_splits} is larger than 1.}
#'     \item{\code{alpha}}{the significance level used.}
#
#'}
#' @export
summary.tsci <- function(object,
                         extended_output = FALSE,
                         ...) {
  stopifnot(inherits(object, "tsci"))


  if (!is.null(object$mult_split_method)) {
    if (object$mult_split_method == "FWER") {
      sd_sel <- rep(".", length(object$Coef_sel))
      sd_all <- rep(".", length(object$Coef_all))
    } else {
      sd_sel <- object$sd_sel
      sd_all <- object$sd_all
    }
  } else {
    sd_sel <- object$sd_sel
    sd_all <- object$sd_all
  }

  coef_df <-
    data.frame(object$Coef_sel[1],
               sd_sel[1],
               object$CI_sel[1, 1],
               object$CI_sel[2, 1],
               object$pval_sel[1])
  colnames(coef_df) <- c("Estimate", "Std_Error", paste(100 * object$alpha/2, "%"), paste(100*(1 - object$alpha/2), "%"), "Pr(>|t|)")
  row.names(coef_df) <- names(object$Coef_sel[1])

  TreatmentModel_df <- data.frame("Estimation_Method" = object$FirstStage_model)

  ViolationSpace_DF <- data.frame("q_comp" = object$q_comp,
                                  "q_cons" = object$q_cons,
                                  "Qmax" = object$Qmax)
  row.names(ViolationSpace_DF) <- names(object$Qmax)

  if (extended_output) {
    coef_all_df <-
      data.frame(object$Coef_all,
                 sd_all,
                 object$CI_all[1, ],
                 object$CI_all[2, ],
                 object$pval_all)
    colnames(coef_all_df) <- c("Estimate", "Std_Error", paste(100 * object$alpha/2, "%"),
                               paste(100*(1 - object$alpha/2), "%"), "Pr(>|t|)")
    row.names(coef_all_df) <- names(object$Coef_all)

    IVStrength_df <- data.frame("IV_Strength" = object$iv_str,
                                "IV_Threshold" = object$iv_thol)
    row.names(IVStrength_df) <- names(object$iv_str)

    output_list <- list(coefficient = coef_df,
                        coefficients_all = coef_all_df,
                        iv_strength = IVStrength_df,
                        viospace_selection = ViolationSpace_DF,
                        treatment_model = TreatmentModel_df,
                        sample_size_A1 = object$n_A1,
                        sample_size_A2 = object$n_A2,
                        n_splits = object$nsplits,
                        mult_split_method = object$mult_split_method,
                        alpha = object$alpha,
                        extended_output = extended_output)
    class(output_list) <- "summary.tsci"
    return(output_list)
  }


  output_list <- list(coefficient = coef_df,
                      invalidity = object$invalidity,
                      viospace_selection = ViolationSpace_DF,
                      treatment_model = TreatmentModel_df,
                      sample_size_A1 = object$n_A1,
                      sample_size_A2 = object$n_A2,
                      n_splits = object$nsplits,
                      mult_split_method = object$mult_split_method,
                      alpha = object$alpha,
                      extended_output = extended_output)
  class(output_list) <- "summary.tsci"
  return(output_list)
}
