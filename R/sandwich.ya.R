#' Calculates standard error of OR estimates of the parameters E(Y1) and E(Y0)
#'
#' This is a support function for \code{\link{ui.causal}} and \code{\link{ui.missing}} and calculates standard error of average of the potential outcomes for the regression imputation estimator.
#' @param X Covariate matrix.
#' @param y Outcome vector.
#' @param z missingness indicator.
#' @param Beta Coefficients from the regression.
#' @param NaivEst Naiv estimates.
#' @param N Total number.
#' @param p Number of covariates outcome regression
#' @param Y1 Binary variable which specifies if standard error corresponds to the estimate of average of the potential outcome under treatment.
#'
#' @importFrom stats binomial coef complete.cases cov get_all_vars glm lm model.matrix pnorm qnorm
#' @export
sandwich.ya <- function(X, y, z, Beta, NaivEst, N, p, Y1 = FALSE) {
  if (Y1) {
    a12 <- apply(X, 2, function(x) {
      mean(-x)
    })
  } else {
    a12 <- apply(X, 2, function(x) {
      mean(x)
    })
  }

  D <- matrix(nrow = (p + 1), ncol = (p + 1), data = 0)
  if (Y1) {
    D <- t(z * X) %*% (z * X) / N
  } else {
    D <- t((1 - z) * X) %*% ((1 - z) * X) / N
  }

  Aninv <- c(1, -a12 %*% solve(D))


  phi <- matrix(nrow = N, ncol = p + 2)
  phi[, 1] <- X %*% Beta - NaivEst

  if (Y1) {
    phi[, 2:(p + 2)] <- apply(X, 2, function(x) {
      z * (y - X %*% Beta) * x
    })
  } else {
    phi[, 2:(p + 2)] <- apply(X, 2, function(x) {
      (1 - z) * (y - X %*% Beta) * x
    })
  }
  # phi[,2:(p+2)]<-ifelse(Y1,apply(X,2,function(x){z*(y-X%*%Beta)*x}),
  #                       apply(X,2,function(x){(1-z)*(y-X%*%Beta)*x}))


  Bn <- matrix(nrow = (2 + p), ncol = (2 + p))
  for (i in 1:(2 + p)) {
    for (j in 1:(2 + p)) {
      Bn[i, j] <- mean(phi[, i] * phi[, j])
    }
  }
  Bn <- Bn / N

  return(Aninv %*% Bn %*% Aninv)
}
