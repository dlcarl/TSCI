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
  if (object$sample_size_A2 == 0) {
    cat(paste("Sample size:", object$sample_size_A1, "\n"))
    cat("No sample splitting was performed.\n")
  } else {
    cat(paste("Sample size A1:", object$sample_size_A1, "\n"))
    cat(paste("Sample size A2:", object$sample_size_A2, "\n"))
    cat(paste("Number of data splits:", object$n_splits, "\n"))
    cat(paste("Aggregation method:", object$mult_split_method, "\n"))
  }

  # print coefficients

  cat("\nCoefficients:\n")
  print(object$coefficients)

  # print treatment model and outcome model statistics

  cat("\nStatistics about the treatment model:\n")
  cat(paste("Estimation method:", object$treatment_model$Estimation.Method, "\n"))
  cat(paste("Residual standard error:", object$treatment_model$Residual.Standard.Error, "\n"))
  cat(paste("R-squared:", object$treatment_model$R.Squared, "\n"))

  cat("\nStatistics about the outcome model:\n")
  print(object$outcome_model)

  cat("\nStatistics about the violation space selection:\n")
  print(object$viospace_selection)

  cat("\nStatistics about the IV strength:\n")
  print(object$iv_strength)

  invisible(object)

}
