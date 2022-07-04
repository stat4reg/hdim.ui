#' Uncertainty intervals for a variable with missingness not at random
#'
#' This function allows you to derive uncertainty intervals for E(Y) with Y being an outcome with missing values not at random. The function uses a regression imputation estimator and a doubly robust estimator. The uncertainty intervals can be used as a sensitivity analysis to unconfoundedness. Note that \code{rho}=0 render the same results as assuming all confounders of Y and missingness indicator variable are observed.
#' @param out.formula Formula for the outcome regression model.
#' @param mis.formula Formula for missingness mechanism (should start with z~). If NULL the same covariates as in the outcome regression will be used.
#' @param data data.frame containing the variables in the formula.
#' @param Y Outcome vector. The arguments Y and X will be only used if out.formula or mis.formula or data matrix is not specified.
#' @param X Covariate matrix. The arguments Y and X will be only used if out.formula or mis.formula or data matrix is not specified.
#' @param rho Pre-specified value of \code{rho}.
#' @param subset Specifies which subset of the covariates should be used in the regression model. It can be "noselection" (default value), "single", "double" or "refit". "Noselection" or "double" is suggested. See details.
#' @param sand Specifies which estimator of the standard errors should be used for OR, see details.
#' @param gridn Number of fixed points within the \code{rho}  intervals for which sigma should be estimated.
#' @param rho.plotrange  An interval larger than \code{rho}. This specifies the range of value (on x axis) where the ui plot will be made using \code{\link{plot.ui}}.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param sigma_correction A character variable which specifies if a corrected estimation of variance of the outcome model is used to find the confounding bias and if so whether we use the old or new correction method. It can take values 'non', 'old' or 'new'.
#'
#' @details In order to visualize the results, you can use \code{\link{plot.ui}}. See the manual for \code{\link{ui}} for details and output.
#'
#' @examples
#' library(MASS)
#' n <- 500
#' rho <- 0.3
#' x <- rnorm(n)
#' s <- 3
#' error <- mvrnorm(n, c(0, 0), matrix(c(1, rho^2, rho^2, s^2), ncol = 2))
#' delta <- c(-0.3, 0.65)
#' zstar <- cbind(1, x) %*% delta + error[, 1]
#' z <- zstar > 0
#' y <- -2 + x - 3 * x^2 + error[, 2]
#'
#'
#' mean(y)
#' y[z==0] <- NA
#' data <- data.frame(y, x)
#'
#' ui <- ui.missing(out.formula =y ~ x + I(x^2), data = data, subset = "double", rho = c(0, 0.4))
#' ui
#' plot(ui)
#' @export

ui.missing <- function(out.formula = NULL,
                       mis.formula = NULL,
                       data = NULL, Y = NULL, X = NULL,
                       rho = c(-0.1, 0.1),
                       subset = "noselection",
                       sand = TRUE, gridn = 21,
                       rho.plotrange = c(-0.5, 0.5), alpha = 0.05,
                       regularization_alpha = 1,
                       sigma_correction = TRUE, ...) {
  treat.formula <- NULL
  if (!is.null(out.formula) & !is.null(data)) {
    names <- all.vars(out.formula)
    Yname <- names[1]
    if (is.null(mis.formula)) {


      char <- as.character(out.formula)
      mis.formula <- as.formula(paste("z~", char[!char %in% c("~", Yname)]))
    }
    treat.formula <- mis.formula
    z <- !is.na(data[, Yname])
    data[!z, Yname] <- 1
    data[, "z"] <- z
  } else {
    z <- !is.na(Y)
    Y[!z] <- 1
  }


  param <- "Y1"
  output <- ui.causal(
    out.formula = out.formula, treat.formula = treat.formula,
    data = data, Y = Y, T = z, X = X,
    rho1 = rho,
    param = param, subset = subset,
    sand = sand, gridn = gridn,
    rho.plotrange = rho.plotrange, alpha = alpha,
    regularization_alpha = regularization_alpha,
    sigma_correction = sigma_correction,
    missing = TRUE
  )


  return(output)
}
