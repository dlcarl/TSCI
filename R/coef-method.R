#' @exportS3Method
coef.tsci <- function(object,
                      all = FALSE,
                      ...) {
  stopifnot(inherits(object, "tsci"))
  if (all) return(c(object$Coef_sel, object$Coef_all))
  return(object$Coef_sel[2])
}
