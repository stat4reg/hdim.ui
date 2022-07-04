#' Correction of OLS sigma for residual of E(Y0|X) needed for derivation of confounding bias
#'
#' This function is a bias correction of the estimator of residual standard deviation of E(Y0|X) under MNAR, for \code{\link{ui.causal}}.
#'
#' @param sigmaOLS Residual sd from outcome regression.
#' @param u Fitted values from propensity score regression.
#' @param gridrho Values of rho.
#' @export

sigmaOLScor0new <- function(sigmaOLS,u, gridrho) {
  
  K = mean(u*lambda0(u))-mean(lambda0(u)^2)
  sqrt(sigmaOLS^2 / (1 + gridrho^2 * K))
  
  
}
