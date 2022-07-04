#' Inverse mills ratio d(x)/(1-p(x)) #' Inverse mills ratio d(x)/p(x) which is the standard normal density function divided by one minus the standard normal cumulative distribution
#'
#' This is a support function for \code{\link{ui.y0}}, \code{\link{ui.y0t1}} and \code{\link{ui.y1t0}} that allows you to calculate the inverse Mills ratio d(x)/(1-p(x)).
#' @param x Vector
#' @importFrom stats pnorm dnorm
#' @export
lambda0 <- function(x) {
  return(dnorm(x) / (1 - pnorm(x)))
}
