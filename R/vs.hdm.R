#' #' Covariate selection using the hdm package
#'
#' This is a support function for \code{\link{ui.causal}} and \code{\link{ui.missing}} and performs variable selection using logistic or linear lasso regression depending on if the output variable is binary or not.
#' @param X Covariate matrix.
#' @param Y outcome vector.
#' @import hdm
#' @export

vs.hdm <- function(Y, X, minscreen = 0) {
  X <- as.matrix(X)
  fit <- ifelse(length(unique(Y)) > 2, rlasso(Y ~ X), rlassologit(Y ~ X))
  return(which((fit[[1]])[-1] != 0))
}
