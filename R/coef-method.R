#' @exportS3Method
coef.tsci <- function(object,
                      all = FALSE,
                      ...) {
  stopifnot(inherits(object, "tsci"))
  if (all) return(c(object$Coef_robust, object$Coef_all))
  return(object$Coef_robust[2])
}
