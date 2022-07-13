#' @exportS3Method
coef.tsci <- function(object,
                      ...) {
  stopifnot(inherits(object, "tsci"))

  betas <- c(object$Coef_robust, object$Coef_all)

  return(betas)
}
