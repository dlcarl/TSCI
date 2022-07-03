#' Title
#'
#' @param object xxx
#' @param ... xxx
#'
#' @return xxx
#' @export
print.summary.tsci <- function(object, ...) {

  # print multi-splitting statistics
  cat("\nStatistics about the data splitting procedure:\n")
  if (object$split_prop == 1) {
    cat("No sample splitting was performed.\n")
  } else {
    cat(paste("Splitting proportion:", round(object$split_prop, 2), "\n"))
    cat(paste("Number of data splits:", object$nsplits, "\n"))
    cat(paste("Aggregation method:", object$mult_split_method, "\n"))
  }

  # print coefficients

  cat("\nCoefficients of selected violation space candidates:\n")
  if (!is.null(object$mult_split_method)) {
    if (object$mult_split_method == "FWER") {
      sd_robust <- t_values <- rep("-", length(object$coef_robust))
    } else {
      t_values <- round(object$coef_robust / object$sd_robust, 3)
      sd_robust <- round(object$sd_robust, 5)
    }
  } else {
    t_values <- round(object$coef_robust / object$sd_robust, 3)
    sd_robust <- round(object$sd_robust, 5)
  }

  coef_df <-
    data.frame(round(object$coef_robust, 5),
               sd_robust,
               t_values,
               round(object$pval_robust, 3))
  colnames(coef_df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(coef_df) <- names(object$coef_robust)
  print(coef_df)

  cat("\nCoefficients of all violation space candidates:\n")
  if (!is.null(object$mult_split_method)) {
    if (object$mult_split_method == "FWER") {
      sd_all <- t_values <- rep("-", length(object$coef_all))
    } else {
      t_values <- round(object$coef_all / object$sd_all, 3)
      sd_all <- round(object$sd_all, 5)
    }
  } else {
    t_values <- round(object$coef_all / object$sd_all, 3)
    sd_all <- round(object$sd_all, 5)
  }
  coef_df <-
    data.frame(round(object$coef_all, 5),
               sd_all,
               t_values,
               round(object$pval_all, 3))
  colnames(coef_df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(coef_df) <- names(object$coef_all)
  print(coef_df)

  # print treatment model and outcome model statistics

  cat("\nStatistics about the treatment model:\n")
  cat(paste("Estimation method:", object$FirstStage_model, "\n"))
  cat(paste("Residual standard error:", round(object$FirstStage_rse, 4), "\n"))
  cat(paste("R-squared:", round(object$FirstStage_Rsquared, 4), "\n"))

  cat("\nStatistics about the outcome model:\n")
  matrix_SecondStage <- rbind(round(object$SecondStage_rse, 4), round(object$SecondStage_Rsquared, 4))
  colnames(matrix_SecondStage) <- names(object$SecondStage_rse)
  rownames(matrix_SecondStage) <- c("Residual standard error", "R-squared")
  print(t(matrix_SecondStage))

  cat("\nStatistics about the violation space selection:\n")
  matrix_viospace <-
    rbind(object$q_comp, object$q_robust, object$Qmax)
  colnames(matrix_viospace) <- names(object$Qmax)
  rownames(matrix_viospace) <- c("q comp", "q robust", "Qmax")
  print(t(matrix_viospace))

  invisible(object)

}
