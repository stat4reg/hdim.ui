#' Print function for intervals
#'
#' This is a support function for \code{\link{print.ui}} which allows you to print an interval (vector of two elements) using parenthesis.
#' @param v Lower and upper bounds.
#' @param digits Number of decimals.
#'
#'
#' @export


interv.p <- function(v, digits = 3) {
  return(paste("(", round(v[1], digits), ", ", round(v[2], digits), ")", sep = ""))
}
