#' Inverse mills ratio d(x)/p(x) which is the standard normal density function divided by the standard normal cumulative distribution
#'
#' This is a support function for \code{\link{ui.y1}}, \code{\link{ui.y1t0}} and \code{\link{ui.y0t1}} that allows you to calculate the inverse Mills ratio d(x)/p(x).
#' @param x Vector
#' @importFrom stats pnorm dnorm
#' @export
lambda1 <- function(x) {
  return(dnorm(x) / pnorm(x))
}
