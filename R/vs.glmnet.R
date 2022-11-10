#' Covariate selection using the glmnet package
#'
#' This is a support function for \code{\link{ui.causal}} and \code{\link{ui.missing}} and performs covariate selection using probit or linear lasso regression depending on if the output variable is binary or not. The lasso regression is preformed using the glmnet package where the regularization parameter is found by cross-validation.
#' @param X Covariate matrix.
#' @param Y outcome vector.
#' @param alpha Parameter that passes to the glmnet function. It is for the elastic net mixing parameter with range [0,1]. alpha=1 is lasso regression (default) and alpha=0 is ridge regression.
#' @import glmnet
#' @export


vs.glmnet <- function(Y, X, alpha, lambda) {
  X <- as.matrix(X)

  if (length(unique(Y)) > 2) {
    fit <- cv.glmnet(x = X, y = as.matrix(Y), family = "gaussian", alpha = alpha)
  } # fit=rlasso(Y~X)
  else {
    fit <- cv.glmnet(
      x = X, y = as.matrix(Y),
      family = binomial(link = "probit"), alpha = alpha
    )
  }

  if(lambda=='lambda1se')
    return(which(coef(fit)[-1] != 0))
  else if(lambda=='lambdamin')
    return(which(coef(fit, s = fit$lambda.min)[-1] != 0))
}
