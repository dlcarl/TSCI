#' Title
#'
#' @param x xxx
#' @param ... xxx
#'
#' @return xxx
#' @export
print.summary.tsci <- function(x, ...) {

  # print multi-splitting statistics
  cat("\nStatistics about the data splitting procedure:\n")
  if (x$sample_size_A2 == 0) {
    cat(paste("Sample size:", x$sample_size_A1, "\n"))
    cat("No sample splitting was performed.\n")
  } else {
    cat(paste("Sample size A1:", x$sample_size_A1, "\n"))
    cat(paste("Sample size A2:", x$sample_size_A2, "\n"))
    cat(paste("Number of data splits:", x$n_splits, "\n"))
    cat(paste("Aggregation method:", x$mult_split_method, "\n"))
  }

  # print coefficients

  cat("\nCoefficients:\n")
  print(x$coefficients)

  # print treatment model and outcome model statistics

  cat("\nStatistics about the treatment model:\n")
  cat(paste("Estimation method:", x$treatment_model$Estimation.Method, "\n"))
  cat(paste("Residual standard error:", x$treatment_model$Residual.Standard.Error, "\n"))
  cat(paste("R-squared:", x$treatment_model$R.Squared, "\n"))

  cat("\nStatistics about the outcome model:\n")
  print(x$outcome_model)

  cat("\nStatistics about the violation space selection:\n")
  print(x$viospace_selection)

  cat("\nStatistics about the IV strength:\n")
  print(x$iv_strength)

  invisible(x)

}
