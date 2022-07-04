#' Correction of OLS sigma for residual of E(Y1|X) needed for derivation of confounding bias
#'
#' This function is a bias correction of the estimator of residual standard deviation of E(Y1|X) under MNAR, for \code{\link{ui.causal}} and \code{\link{ui.missing}}.
#'
#' @param sigmaOLS Residual sd from outcome regression.
#' @param u Fitted values from propensity score regression.
#' @param gridrho Values of rho.
#' @export
# sigmaOLScor1(Xy.design[t==1,],sigma1hatOLS,n1,p,u1,gridrho1)
sigmaOLScor1new <- function( sigmaOLS, u, gridrho) {

  K = - mean(u*lambda1(u))-mean(lambda1(u)^2)
  
  sqrt(sigmaOLS^2 / (1 + gridrho^2 * K))
}
