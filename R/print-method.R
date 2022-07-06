#' Title xxx
#'
#' @param x xxx
#' @param ... xxx
#'
#' @return xxx
#' @exportS3Method
print.tsci <- function(x, ...) {
  stopifnot(inherits(x, "tsci"))

  print.summary.tsci(summary.tsci(object = x, ...), ...)
}
