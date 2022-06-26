#' Title
#'
#' @param list_outputs xxx
#' @param ind_start xxx
#'
#' @return xxx
#' @noRd
#'
check_output <- function(list_outputs, ind_start) {
  n <- length(list_outputs)
  ind_na <- sapply(list_outputs, FUN = function(x) any(is.na(x$value)))
  prop_na <- sum(ind_na) / n
  errors <- lapply(list_outputs, FUN = function(x) x$error)
  messages <- lapply(list_outputs, FUN = function(x) x$message)
  warnings <- lapply(list_outputs, FUN = function(x) x$warning)
  error_string <- NULL
  for (i in seq_len(n)) {
    if(!is.null(errors[[i]]) | !is.null(messages[[i]]) | !is.null(warnings[[i]])){
      error_string <- paste(error_string, paste("split", i + ind_start - 1),
                            "The following errors occured:",
                            errors[[i]],
                            "The following warnings were raised:",
                            warnings[[i]],
                            "The following messages were raised:",
                            messages[[i]],
                            sep = "\n")
    }
  }

  return(list(prop_na  = prop_na,
              ind_na = ind_na,
              error_string = error_string))
}
