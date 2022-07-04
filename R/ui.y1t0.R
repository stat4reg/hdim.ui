#' Uncertainty intervals for E(Y1|T=0)
#'
#' This is a support function for \code{\link{ui.causal}} and derives the uncertainty intervals for E(Y1|T=0).
#' @param out.formula Formula for the outcome regression model
#' @param y.data data.frame containing variables needed for fitting the outcome model
#' @param gamma Vector of estimated coefficients from the fitted model for treatment variable
#' @param t.data data.frame containing variables needed for fitting the propensity score model
#' @param rho1 Pre-specified interval for \code{rho1}.
#' @param sand Specifies which estimator of the standard errors should be used for OR, see details in \code{\link{ui.causal}}.
#' @param gridn Number of fixed points within the \code{rho1} interval for which sigma1 should be estimated.
#' @param rho.plotrange an interval larger than \code{rho1} for the plot using \code{\link{plot.ui}}.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param sigma_correction A character variable which specifies if a corrected estimation of variance of the outcome model is used to find the confounding bias and if so whether we use the old or new correction method. It can take values 'non', 'old' or 'new'.
#'
#' @importFrom stats binomial coef complete.cases cov get_all_vars glm lm model.matrix pnorm qnorm
#' @export
ui.y1t0 <- function(out.formula, y.data,
                    gamma, t.data,
                    rho1 = c(-0.3, 0.3),
                    sand = TRUE, gridn = 21,
                    rho.plotrange = c(-0.5, 0.5), alpha = 0.05, sigma_correction, ...) {
  plot <- TRUE
  output <- list()
  output$plot <- list()
  output$plot$plot <- plot

  t <- t.data[, 1]
  y <- y.data[, 1]
  y1 <- y[t == 1]


  n0 <- sum(t == 0)
  n1 <- sum(t == 1)
  N <- n0 + n1
  p <- ncol(y.data) - 1

  out.model1 <- lm(out.formula, data = y.data[t == 1, ])
  output$out.model1 <- out.model1

  Xy.design <- model.matrix(lm(out.formula, data = y.data))
  Xt.design <- model.matrix(as.formula(paste("~", paste(names(t.data)[-1], collapse = "+"))), data = t.data)


  # Estimating BetaOLS treated
  BetaOLSy1 <- coef(summary(out.model1))[, 1]
  varBetaOLSy1 <- diag((coef(summary(out.model1))[, 2])^2)
  sigma1hatOLS <- summary(out.model1)$sigma

  u <- (Xt.design %*% gamma)
  u1 <- (Xt.design[t == 1, ] %*% gamma)
  OLSlambda1 <- solve(t(Xy.design[t == 1, ]) %*% Xy.design[t == 1, ]) %*% t(Xy.design[t == 1, ]) %*% lambda1(u1)
  l1 <- mean(lambda1(u))

  u0 <- (Xt.design[t == 0, , drop = FALSE] %*% gamma)
  l0t0 <- mean(lambda0(u0))


  t1 <- gridrho.f(rho1, gridn, rho.plotrange, plot)
  gridrho1 <- t1[[1]]
  nui1 <- t1[[2]]
  output$plot$nui1 <- nui1



  if (sigma_correction=='non') { output$plot$sigma1 <- sigma1hatOLS
  } else if (sigma_correction=='old') {output$plot$sigma1 <- sigmaOLScor1(Xy.design[t == 1, ], sigma1hatOLS, n1, p, u1, gridrho1)
  }else if (sigma_correction=='new'){output$plot$sigma1 <- sigmaOLScor1new(sigma1hatOLS,u1, gridrho1)
  }


  output$sigma1 <- output$plot$sigma1[nui1]
  output$rho1 <- rho1


  output$plot$gridrho1 <- gridrho1
  output$plot$gridrho <- gridrho1

  output$gridrho1 <- gridrho1[nui1]


  output$DR$conf.bias <- vector(length = 2)
  names(output$DR$conf.bias) <- c("Min", "Max")
  output$DR$IdentInt <- output$DR$conf.bias
  output$DR$ui <- output$DR$conf.bias
  output$DR$naivci <- output$DR$conf.bias
  output$OR <- output$DR

  ##################### NAIV ESTIMATES #####################
  phat <- pnorm(Xt.design %*% gamma)
  phat1 <- phat[t == 1]

  output$OR$naiv <- mean(Xy.design[t == 0, ] %*% BetaOLSy1)
  output$DR$naiv <- output$OR$naiv + sum((y1 - Xy.design[t == 1, ] %*% BetaOLSy1) / (phat1)) / n0


  ##################### CONFOUNDING BIAS #####################
  exgt0 <- rep(1, n0) %*% Xy.design[t == 0, ] / n0


  # Counfounding bias regression imputation ACT
  f <- gridrho1 * output$plot$sigma1 * c(l0t0 + exgt0 %*% OLSlambda1)
  names(f) <- as.character(round(gridrho1, 4))
  output$plot$OR$conf.bias.grid <- f
  output$OR$conf.bias.grid <- f[nui1]
  output$OR$conf.bias <- c(min(f[nui1]), max(f[nui1]))
  # Total bias DR ACT
  f <- gridrho1 * output$plot$sigma1 * l1 * N / n0
  names(f) <- as.character(round(gridrho1, 4))
  output$plot$DR$conf.bias.grid <- f
  output$DR$conf.bias.grid <- f[nui1]
  output$DR$conf.bias <- c(min(f[nui1]), max(f[nui1]))

  # Coef, NAIV ESTIMATES - BIAS (for each value of rho)
  # Want the same rownames and colnames as conf.bias.grid
  output$plot$DR$coef <- output$plot$DR$conf.bias.grid
  output$plot$DR$coef <- c(output$DR$naiv) - output$plot$DR$conf.bias.grid
  output$plot$OR$coef <- output$plot$OR$conf.bias.grid
  output$plot$OR$coef <- c(output$OR$naiv) - output$plot$OR$conf.bias.grid
  output$DR$coef <- output$DR$conf.bias.grid
  output$DR$coef <- c(output$DR$naiv) - output$DR$conf.bias.grid
  output$OR$coef <- output$OR$conf.bias.grid
  output$OR$coef <- c(output$OR$naiv) - output$OR$conf.bias.grid

  # Identification interval, NAIV ESTIMATES - BIAS (interval)
  output$OR$IdentInt <- rep(output$OR$naiv, 2) - output$OR$conf.bias[2:1]
  output$DR$IdentInt <- rep(output$DR$naiv, 2) - output$DR$conf.bias[2:1]


  # Standard errors
  if (sand == TRUE) {
    output$OR$se <- sqrt(sandwich.yatb(Xy.design, y, t, BetaOLSy1, output$OR$naiv, N, p, Y1 = TRUE))
  } else {
    output$OR$se <- sqrt(((BetaOLSy1) %*% cov(Xy.design[t == 0, ]) %*% (BetaOLSy1)) / n0)
  }

  I <- ((1 - t) * Xy.design %*% BetaOLSy1 + (t) * (1 - phat) * (y - Xy.design %*% BetaOLSy1) / (phat)) / n0 * N - output$DR$naiv
  output$DR$se <- sqrt(sum(I^2)) / N


  zalpha <- qnorm(1 - alpha / 2)
  # Confidence intervals

  dim <- c(length(output$plot$OR$conf.bias.grid), 2)
  dnames <- list(names(output$plot$OR$conf.bias.grid), c("lower", "upper"))

  output$plot$OR$ci <- array(c(c(output$OR$naiv - zalpha * output$OR$se) - output$plot$OR$conf.bias.grid, c(output$OR$naiv + zalpha * output$OR$se) -
    output$plot$OR$conf.bias.grid), dim = dim, dimnames = dnames)
  output$plot$DR$ci <- array(c(c(output$DR$naiv - zalpha * output$DR$se) - output$plot$DR$conf.bias.grid, c(output$DR$naiv + zalpha * output$DR$se) -
    output$plot$DR$conf.bias.grid), dim = dim, dimnames = dnames)

  dim <- c(length(output$OR$conf.bias.grid), 2)
  dnames <- list(names(output$OR$conf.bias.grid), c("lower", "upper"))

  output$OR$ci <- array(c(c(output$OR$naiv - zalpha * output$OR$se) - output$OR$conf.bias.grid, c(output$OR$naiv + zalpha * output$OR$se) -
    output$OR$conf.bias.grid), dim = dim, dimnames = dnames)
  output$DR$ci <- array(c(c(output$DR$naiv - zalpha * output$DR$se) - output$DR$conf.bias.grid, c(output$DR$naiv + zalpha * output$DR$se) -
    output$DR$conf.bias.grid), dim = dim, dimnames = dnames)

  output$OR$ui <- c(min(output$OR$ci), max(output$OR$ci))
  output$DR$ui <- c(min(output$DR$ci), max(output$DR$ci))
  output$OR$naivci <- output$OR$naiv + c(-zalpha * output$OR$se, zalpha * output$OR$se)
  output$DR$naivci <- output$DR$naiv + c(-zalpha * output$DR$se, zalpha * output$DR$se)
  output$plot$gridn <- gridn

  names(output$OR$ui) <- c("lower", "upper")
  names(output$OR$naivci) <- c("lower", "upper")
  names(output$OR$conf.bias) <- c("lower", "upper")
  names(output$OR$IdentInt) <- c("lower", "upper")
  names(output$DR$ui) <- c("lower", "upper")
  names(output$DR$naivci) <- c("lower", "upper")
  names(output$DR$conf.bias) <- c("lower", "upper")
  names(output$DR$IdentInt) <- c("lower", "upper")


  output$call <- match.call()
  class(output) <- "ui"

  return(output)
}
