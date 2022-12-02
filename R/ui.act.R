#' Uncertainty intervals for average causal effect on treated
#'
#' This is a support function for \code{\link{ui.causal}} and derives uncertainty intervals for average causal effect on treated.
#' @param out.formula Formula for the outcome regression model
#' @param treat.formula Formula for the propensity score model (regression model for treatment assignment).
#' @param data data.frame containing the variables in the formulas including dummy variables.
#' @param gamma Vector of estimated coefficients from the fitted model for treatment variable
#' @param rho0 Pre-specified interval for \code{rho0}.
#' @param sand Specifies which estimator of the standard errors should be used for OR, see details in \code{\link{ui.causal}}.
#' @param gridn Number of fixed points within the \code{rho0} interval for which sigma1 should be estimated.
#' @param rho.plotrange an interval larger than \code{rho0} for the plot using \code{\link{plot.ui}}.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param sigma_correction A character variable which specifies if a corrected estimation of variance of the outcome model is used to find the confounding bias and if so whether we use the old or new correction method. It can take values 'non', 'old' or 'new'.
#'
#' @importFrom stats binomial coef complete.cases cov get_all_vars glm lm model.matrix pnorm qnorm
#' @export
ui.act <- function(out.formula, treat.formula, data, gamma,
                   rho0 = c(-0.3, 0.3),
                   sand = TRUE, gridn = 21,
                   rho.plotrange = c(-0.5, 0.5), alpha = 0.05, sigma_correction, ...) {
  output0 <- ui.y0t1(out.formula, treat.formula, data, gamma,
    rho0,
    sand, gridn,
    rho.plotrange, alpha, sigma_correction
  )

  t <- t.data[, 1]
  y <- y.data[, 1]
  y1 <- y[t == 1]

  n1 <- sum(t == 1)


  plot <- TRUE
  output <- list()
  output$plot <- list()
  output$plot$plot <- plot

  output$out.model0 <- output0$out.model0

  output$plot$nui0 <- output0$plot$nui0
  output$plot$sigma0 <- output0$plot$sigma0
  output$sigma0 <- output$sigma0
  output$rho0 <- rho0

  output$plot$gridrho <- output0$plot$gridrho
  output$plot$gridrho0 <- output0$plot$gridrho0
  output$gridrho0 <- output0$gridrho0


  output$DR$conf.bias <- vector(length = 2)
  names(output$DR$conf.bias) <- c("Min", "Max")
  output$DR$IdentInt <- output$DR$conf.bias
  output$DR$ui <- output$DR$conf.bias
  output$DR$naivci <- output$DR$conf.bias
  output$OR <- output$DR

  ##################### NAIV ESTIMATES #####################
  output$OR$naiv <- mean(y1) - output0$OR$naiv
  output$DR$naiv <- mean(y1) - output0$DR$naiv
  #* important: changed -bias to +bias (because of - in -E(Y0T1))
  output$plot$OR$conf.bias.grid <- -output0$plot$OR$conf.bias.grid
  output$plot$DR$conf.bias.grid <- -output0$plot$DR$conf.bias.grid
  output$OR$conf.bias.grid <- -output0$OR$conf.bias.grid
  output$DR$conf.bias.grid <- -output0$DR$conf.bias.grid
  output$OR$conf.bias <- -output0$OR$conf.bias
  output$DR$conf.bias <- -output0$DR$conf.bias



  output$plot$OR$coef <- mean(y1) - output0$plot$OR$coef
  output$plot$DR$coef <- mean(y1) - output0$plot$DR$coef

  output$OR$coef <- mean(y1) - output0$OR$coef
  output$DR$coef <- mean(y1) - output0$DR$coef

  output$OR$se <- sqrt(sd(y1)^2 / n1 + output0$OR$se^2)
  output$DR$se <- sqrt(sd(y1)^2 / n1 + output0$DR$se^2)

  output$OR$IdentInt <- rev(mean(y1) - output0$OR$IdentInt)
  output$DR$IdentInt <- rev(mean(y1) - output0$DR$IdentInt)



  # output$plot$OR$ci[,1]  <-mean(y1)-output$plot$OR$ci[,2]
  # output$plot$OR$ci[,2]  <-mean(y1)-output$plot$OR$ci[,1]
  # output$plot$DR$ci[,1]  <-mean(y1)-output$plot$DR$ci[,2]
  # output$plot$DR$ci[,2]  <-mean(y1)-output$plot$DR$ci[,1]
  #
  # output$OR$ui <-rev(mean(y1)- output$OR$ui)
  # output$DR$ui <-rev(mean(y1)- output$DR$ui)
  #
  # output$OR$naivci<-rev(mean(y1)-  output$OR$naivci)
  # output$DR$naivci <-rev(mean(y1)- output$DR$naivci)



  zalpha <- qnorm(1 - alpha / 2) # Confidence intervals needs to be calculated taking into account variability in E(y1)


  dim <- c(length(output$plot$OR$conf.bias.grid), 2)
  dnames <- list(names(output$plot$OR$conf.bias.grid), c("lower", "upper"))

  output$plot$OR$ci <- array(c(
    c(output$OR$naiv - zalpha * output$OR$se) - output$plot$OR$conf.bias.grid,
    c(output$OR$naiv + zalpha * output$OR$se) - output$plot$OR$conf.bias.grid
  ), dim = dim, dimnames = dnames)
  output$plot$DR$ci <- array(c(
    c(output$DR$naiv - zalpha * output$DR$se) - output$plot$DR$conf.bias.grid,
    c(output$DR$naiv + zalpha * output$DR$se) - output$plot$DR$conf.bias.grid
  ), dim = dim, dimnames = dnames)


  dim <- c(length(output$OR$conf.bias.grid), 2)
  dnames <- list(names(output$OR$conf.bias.grid), c("lower", "upper"))

  output$OR$ci <- array(c(
    c(output$OR$naiv - zalpha * output$OR$se) - output$OR$conf.bias.grid,
    c(output$OR$naiv + zalpha * output$OR$se) - output$OR$conf.bias.grid
  ), dim = dim, dimnames = dnames)
  output$DR$ci <- array(c(
    c(output$DR$naiv - zalpha * output$DR$se) - output$DR$conf.bias.grid,
    c(output$DR$naiv + zalpha * output$DR$se) - output$DR$conf.bias.grid
  ), dim = dim, dimnames = dnames)

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
