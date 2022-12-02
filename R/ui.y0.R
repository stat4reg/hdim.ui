#' Uncertainty intervals for E(Y0)
#'
#' This is a support function for \code{\link{ui.causal}} and derives the uncertainty intervals for E(Y0)
#' @param out.formula Formula for the outcome regression model
#' @param treat.formula Formula for the propensity score model (regression model for treatment assignment).
#' @param data data.frame containing the variables in the formulas including dummy variables.
#' @param gamma Vector of estimated coefficients from the fitted model for treatment variable
#' @param rho0 Pre-specified interval for \code{rho0}.
#' @param sand Specifies which estimator of the standard errors should be used for OR, see details in \code{\link{ui.causal}}.
#' @param gridn Number of fixed points within the \code{rho0} interval for which sigma0 should be estimated.
#' @param rho.plotrange an interval larger than \code{rho0} for the plot using \code{\link{plot.ui}}.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param sigma_correction A character variable which specifies if a corrected estimation of variance of the outcome model is used to find the confounding bias and if so whether we use the old or new correction method. It can take values 'non', 'old' or 'new'.
#'
#' @importFrom stats binomial coef complete.cases cov get_all_vars glm lm model.matrix pnorm qnorm
#' @export
ui.y0 <- function(out.formula, treat.formula, data, gamma,
                  rho0 = c(-0.3, 0.3),
                  sand = TRUE, gridn = 21,
                  rho.plotrange = c(-0.5, 0.5), alpha = 0.05, sigma_correction, ...) {
  plot <- TRUE
  output <- list()
  output$plot <- list()
  output$plot$plot <- plot

  out.names <- all.vars(out.formula)
  Yname <- out.names[1]
  treat.names <- all.vars(treat.formula)
  Tname <- treat.names[1]

  t <- data[, Tname]
  y <- data[, Yname]
  y0 <- y[t == 0]


  n0 <- sum(t == 0)
  n1 <- sum(t == 1)
  N <- n0 + n1
  p <- length(attr(terms(out.formula), "term.labels"))



  out.model0 <- lm(out.formula, data = data[t == 0, ])
  output$out.model0 <- out.model0

  Xy.design <- model.matrix(out.formula, data = data)
  Xt.design <- model.matrix(treat.formula, data = data)

  # Estimating BetaOLS untreated
  BetaOLSy0 <- coef(summary(out.model0))[, 1]
  varBetaOLSy0 <- diag((coef(summary(out.model0))[, 2])^2)
  sigma0hatOLS <- summary(out.model0)$sigma

  u <- (Xt.design %*% gamma)
  u0 <- (Xt.design[t == 0, , drop = FALSE] %*% gamma)
  OLSlambda0 <- solve(t(Xy.design[t == 0, ]) %*% Xy.design[t == 0, ]) %*% t(Xy.design[t == 0, ]) %*% lambda0(u0)
  l0 <- mean(lambda0(u))


  t0 <- gridrho.f(rho0, gridn, rho.plotrange, plot)
  gridrho0 <- t0[[1]]
  nui0 <- t0[[2]]
  output$plot$nui0 <- nui0

  xb=Xy.design %*% BetaOLSy0
  xb0=Xy.design[t == 0, ] %*% BetaOLSy0

  if (sigma_correction=='non') {
    output$plot$sigma0 <- sigma0hatOLS
  } else if (sigma_correction=='old') {
    output$plot$sigma0 <- sigmaOLScor0(Xy.design[t == 0, ], sigma0hatOLS, n0, p, u0, gridrho0)
  }else if (sigma_correction=='new'){
    output$plot$sigma0 <- sigmaOLScor0new(sigma0hatOLS, u0, gridrho0)
  }
  output$sigma0 <- output$plot$sigma0[nui0]
  output$rho0 <- rho0
  # AA we now need both gridrho is used when plotting and the other is needed when estimating ate
  output$plot$gridrho0 <- gridrho0
  output$plot$gridrho <- gridrho0

  output$gridrho0 <- gridrho0[nui0]


  output$DR$conf.bias <- vector(length = 2)
  names(output$DR$conf.bias) <- c("Min", "Max")
  output$DR$IdentInt <- output$DR$conf.bias
  output$DR$ui <- output$DR$conf.bias
  output$DR$naivci <- output$DR$conf.bias
  output$OR <- output$DR

  ##################### NAIV ESTIMATES #####################
  phat <- pnorm(u)
  phat0 <- phat[t == 0]

  output$OR$naiv <- (sum(Xy.design[t == 1, ] %*% BetaOLSy0) + sum(xb0)) / N
  output$DR$naiv <- output$OR$naiv + (sum((y0 - xb0) / (1 - phat0))) / N

  ##################### CONFOUNDING BIAS #####################
  exgt0 <- rep(1, n0) %*% Xy.design[t == 0, ] / n0

  # Counfounding bias regression imputation ACE
  f2 <- -gridrho0 * output$plot$sigma0 * c(rep(1, N) %*% Xy.design %*% OLSlambda0) / N
  output$plot$OR$conf.bias.grid <- f2
  names(output$plot$OR$conf.bias.grid) <- as.character(round(gridrho0, 4))
  output$OR$conf.bias.grid <- output$plot$OR$conf.bias.grid[nui0]
  output$OR$conf.bias <- c(min(output$OR$conf.bias.grid), max(output$OR$conf.bias.grid))

  # Total bias DR ACE
  f2 <- -gridrho0 * output$plot$sigma0 * l0
  output$plot$DR$conf.bias.grid <- f2
  names(output$plot$DR$conf.bias.grid) <- names(output$plot$OR$conf.bias.grid)
  output$DR$conf.bias.grid <- output$plot$DR$conf.bias.grid[nui0]
  output$DR$conf.bias <- c(min(output$DR$conf.bias.grid), max(output$DR$conf.bias.grid))

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



  if (sand == TRUE) {
    output$OR$se <- sqrt(sandwich.ya(Xy.design, y, t, BetaOLSy0, output$OR$naiv, N, p, Y1 = FALSE))
  } else {
    output$OR$se <- sqrt(rep(1, N) %*% Xy.design %*% (
      # sigma1hatOLS^2*solve(t(Xy.design[t==1,])%*%Xy.design[t==1,])+
      sigma0hatOLS^2 * solve(t(Xy.design[t == 0, ]) %*% Xy.design[t == 0, ])) %*% t(Xy.design) %*% rep(1, N) / N^2 +
      (BetaOLSy0) %*% cov(Xy.design) %*% (BetaOLSy0) / N)
  }
  I <- (xb + (1 - t) * ((y - xb) / (1 - phat))) - output$DR$naiv
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
