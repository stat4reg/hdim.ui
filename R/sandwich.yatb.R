#' Calculates standard error of OR estimates of the parameters E(Y1|T=0) and E(Y0|T=1)
#'
#' This is a support function for \code{\link{ui.causal}} and calculates standard error of the parameters E(Y1|T=0) and E(Y0|T=1) for the regression imputation estimator.
#' @param X Covariate matrix.
#' @param y Outcome vector.
#' @param z missingness indicator.
#' @param Beta Coefficients from the regression.
#' @param NaivEst Naiv estimates.
#' @param N Total number.
#' @param p Number of covariates outcome regression
#' @param Y1 Binary variable which specifies if standard error corresponds to the estimate of conditional average of the potential outcome under treatment.
#'
#' @importFrom stats binomial coef complete.cases cov get_all_vars glm lm model.matrix pnorm qnorm
#' @export
sandwich.yatb <- function(X, y, z, Beta, NaivEst, N, p, Y1 = FALSE) {
  n <- ifelse(Y1, sum(z == 0), sum(z == 1))

  if (Y1) {
    a12 <- apply(X, 2, function(x) {
      mean(-x * (1 - z))
    })
  } else {
    a12 <- apply(X, 2, function(x) {
      mean(z * x)
    })
  }

  D <- matrix(nrow = (p + 1), ncol = (p + 1), data = 0)

  if (Y1) {
    D <- t(z * X) %*% (z * X) / N
  } else {
    D <- t((1 - z) * X) %*% ((1 - z) * X) / N
  }

  Aninv <- N / n * c(1, -a12 %*% solve(D))


  phi <- matrix(nrow = N, ncol = p + 2)
  if (Y1) {
    phi[, 1] <- (1 - z) * (X %*% Beta - NaivEst)
  } else {
    phi[, 1] <- z * (X %*% Beta - NaivEst)
  }
  if (Y1) {
    phi[, 2:(p + 2)] <- apply(X, 2, function(x) {
      z * (y - X %*% Beta) * x
    })
  } else {
    phi[, 2:(p + 2)] <- apply(X, 2, function(x) {
      (1 - z) * (y - X %*% Beta) * x
    })
  }


  Bn <- matrix(nrow = (p + 2), ncol = (p + 2))
  for (i in 1:(p + 2)) {
    for (j in 1:(p + 2)) {
      Bn[i, j] <- mean(phi[, i] * phi[, j])
    }
  }
  Bn <- Bn / N
  Aninv[1] * Bn[1, 1]
  return(Aninv %*% Bn %*% Aninv)
}
