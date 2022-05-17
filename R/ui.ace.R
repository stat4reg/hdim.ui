#' Uncertainty intervals for average causal effect
#'
#' This is a support function for \code{\link{ui.causal}} and derives uncertainty intervals for average causal effect.
#' @param out.formula Formula for the outcome regression models
#' @param y.data data.frame containing variables needed for fitting the outcome models
#' @param gamma Vector of estimated coefficients from the fitted model for treatment variable
#' @param t.data data.frame containing variables needed for fitting the propensity score model
#' @param rho0 Pre-specified interval for \code{rho0}.
#' @param rho1 Pre-specified interval for \code{rho1}.
#' @param sand Specifies which estimator of the standard errors should be used for OR, see details in \code{\link{ui.causal}}.
#' @param gridn Number of fixed points within the \code{rho0} and \code{rho0} intervals for which sigma0 and sigma1 should be estimated.
#' @param rho.plotrange an interval larger than \code{rho0} and \code{rho1} for the plot using \code{\link{plot.ui}}.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param sigma_correction A logical variable which specifies if a corrected estimation of variance of the outcome model is used to find the confounding bias. Default value is TRUE.
#'
#' @importFrom stats binomial coef complete.cases cov get_all_vars glm lm model.matrix pnorm qnorm
#' @export
#'
#'
ui.ace <- function(out.formula, y.data,
                   gamma, t.data,
                   rho0 = c(-0.3, 0.3), rho1 = c(-0.3, 0.3),
                   sand = TRUE, gridn = 21,
                   rho.plotrange = c(-0.5, 0.5), alpha = 0.05, sigma_correction, ...) {
  plot <- TRUE
  if (plot == TRUE) {
    if (sum(rho0 != rho1) > 0) {
      # warning('In order to plot results rho0 and rho1 must be equal intervals, plot changed to FALSE.')
      plot <- FALSE
    }
  }


  output0 <- ui.y0(
    out.formula, y.data,
    gamma, t.data,
    rho0,
    sand, gridn,
    rho.plotrange, alpha, sigma_correction
  )

  output1 <- ui.y1(
    out.formula, y.data,
    gamma, t.data,
    rho1,
    sand, gridn,
    rho.plotrange, alpha, sigma_correction
  )


  output <- list()
  output$plot <- list()
  output$plot$plot <- plot



  output$plot$nui0 <- output0$plot$nui0
  output$plot$sigma0 <- output0$plot$sigma0
  output$sigma0 <- output$sigma0
  output$rho0 <- rho0
  # AA
  # output$plot$gridrho0<-gridrho0
  output$plot$gridrho <- output0$plot$gridrho # or output1, doesnt matter
  output$plot$gridrho0 <- output0$plot$gridrho0
  output$gridrho0 <- output0$gridrho0

  output$plot$nui1 <- output1$plot$nui1
  output$plot$sigma1 <- output1$plot$sigma1
  output$sigma1 <- output1$sigma1
  output$rho1 <- rho1
  # AA
  # output$plot$gridrho1<-gridrho1
  output$plot$gridrho1 <- output1$plot$gridrho1
  output$gridrho1 <- output1$gridrho1




  output$DR$conf.bias <- vector(length = 2)
  names(output$DR$conf.bias) <- c("Min", "Max")
  output$DR$IdentInt <- output$DR$conf.bias
  output$DR$ui <- output$DR$conf.bias
  output$DR$naivci <- output$DR$conf.bias
  output$OR <- output$DR

  ##################### NAIV ESTIMATES #####################
  output$OR$naiv <- output1$OR$naiv - output0$OR$naiv
  output$DR$naiv <- output1$DR$naiv - output0$DR$naiv
  ##################### CONFOUNDING BIAS #####################
  # Counfounding bias regression imputation ACE
  f1 <- output1$plot$OR$conf.bias.grid
  f2 <- -output0$plot$OR$conf.bias.grid
  output$plot$OR$conf.bias.grid <- matrix(rep(f1, length(f2)), nrow = length(f1)) + t(matrix(rep(f2, length(f1)), nrow = length(f2)))
  rownames(output$plot$OR$conf.bias.grid) <- as.character(round(output$plot$gridrho1, 4))
  colnames(output$plot$OR$conf.bias.grid) <- as.character(round(output$plot$gridrho0, 4))

  # output$OR$conf.bias.grid<-output$plot$OR$conf.bias.grid[output1$plot$nui1,output0$plot$nui0]
  # new, to prevent error when rho has one value only
  f1 <- output1$OR$conf.bias.grid
  f2 <- -output0$OR$conf.bias.grid
  output$OR$conf.bias.grid <- matrix(rep(f1, length(f2)), nrow = length(f1)) + t(matrix(rep(f2, length(f1)), nrow = length(f2)))
  rownames(output$OR$conf.bias.grid) <- as.character(round(output$gridrho1, 4))
  colnames(output$OR$conf.bias.grid) <- as.character(round(output$gridrho0, 4))


  output$OR$conf.bias <- c(min(output$OR$conf.bias.grid), max(output$OR$conf.bias.grid))


  # Total bias DR ACE
  f1 <- output1$plot$DR$conf.bias.grid
  f2 <- -output0$plot$DR$conf.bias.grid
  output$plot$DR$conf.bias.grid <- matrix(rep(f1, length(f2)), nrow = length(f1)) + t(matrix(rep(f2, length(f1)), nrow = length(f2)))
  dimnames(output$plot$DR$conf.bias.grid) <- dimnames(output$plot$OR$conf.bias.grid)

  # output$DR$conf.bias.grid<-output$plot$DR$conf.bias.grid[output1$plot$nui1,output0$plot$nui0]
  # new, to prevent error when rho has one value only
  f1 <- output1$DR$conf.bias.grid
  f2 <- -output0$DR$conf.bias.grid
  output$DR$conf.bias.grid <- matrix(rep(f1, length(f2)), nrow = length(f1)) + t(matrix(rep(f2, length(f1)), nrow = length(f2)))
  dimnames(output$DR$conf.bias.grid) <- dimnames(output$OR$conf.bias.grid)


  output$DR$conf.bias <- c(min(output$DR$conf.bias.grid), max(output$DR$conf.bias.grid))

  ## AA
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
  output$OR$se <- sqrt(output1$OR$se^2 + output0$OR$se^2)
  output$DR$se <- sqrt(output1$DR$se^2 + output0$DR$se^2)


  zalpha <- qnorm(1 - alpha / 2)
  # Confidence intervals

  dim <- c(dim(output$plot$OR$conf.bias.grid), 2)
  dnames <- dimnames(output$plot$OR$conf.bias.grid)
  dnames[[3]] <- c("lower", "upper")

  output$plot$OR$ci <- array(c(
    c(output$OR$naiv - zalpha * output$OR$se) - output$plot$OR$conf.bias.grid,
    c(output$OR$naiv + zalpha * output$OR$se) - output$plot$OR$conf.bias.grid
  ), dim = dim, dimnames = dnames)
  output$plot$DR$ci <- array(c(
    c(output$DR$naiv - zalpha * output$DR$se) - output$plot$DR$conf.bias.grid,
    c(output$DR$naiv + zalpha * output$DR$se) - output$plot$DR$conf.bias.grid
  ), dim = dim, dimnames = dnames)


  dim <- c(dim(output$OR$conf.bias.grid), 2)
  dnames <- dimnames(output$OR$conf.bias.grid)
  dnames[[3]] <- c("lower", "upper")

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
