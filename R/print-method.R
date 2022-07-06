#' Title xxx
#'
#' @param x xxx
#' @param ... xxx
#'
#' @return xxx
#' @export
#'
print.tsci <- function(x, ...) {
  stopifnot(inherits(x, "tsci"))

  print.summary.tsci(summary.tsci(object = x, ...), ...)
}
