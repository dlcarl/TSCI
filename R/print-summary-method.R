#' @exportS3Method
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
  cat("\nTreatment effect estimate of selected violation space candidate:\n")
  coefficients_df <- data.frame(lapply(x$coefficients, FUN = function(y) if (is.numeric(y)) {round(y, 5)}  else {y}))
  rownames(coefficients_df) <- rownames(x$coefficients)
  colnames(coefficients_df) <- colnames(x$coefficients)
  print(coefficients_df)

  if (x$extended_output) {
    cat("\nTreatment effect estimates of all violation space candidates:\n")
    coefficients_all_df <- data.frame(lapply(x$coefficients_all, FUN = function(y) if (is.numeric(y)) {round(y, 5)}  else {y}))
    rownames(coefficients_all_df) <- rownames(x$coefficients_all)
    colnames(coefficients_all_df) <- colnames(x$coefficients_all)
    print(coefficients_all_df)
  }

  # print treatment model and outcome model statistics

  cat("\nStatistics about the treatment model:\n")
  cat(paste("Estimation method:", x$treatment_model$Estimation_Method, "\n"))
  # cat(paste("Residual standard error:", round(x$treatment_model$Residual_Standard_Error, 4), "\n"))
  # cat(paste("R-squared:", round(x$treatment_model$R_Squared, 4), "\n"))

  # cat("\nStatistics about the outcome model:\n")
  # print(round(x$outcome_model, 4))

  cat("\nStatistics about the violation space selection:\n")
  print(x$viospace_selection)

  if (x$extended_output) {
    cat("\nStatistics about the IV strength:\n")
    print(round(x$iv_strength, 2))
  }

  invisible(x)

}
