#' Uncertainty intervals for causal parameters including average causal effect and average causal effect among treated
#'
#' This function allows you to derive uncertainty intervals for causal parameters E(Y1), E(Y0), E(Y1|T=0), E(Y0|T=1) in an observational study with Yt being the potential outcome that would have been observed under treatment T=t. The function uses a regression imputation estimator (OR) and a doubly robust estimator (DR). The uncertainty intervals can be used as a sensitivity analysis to unconfoundedness. Note that \code{rho0}=0 and \code{rho1}=0 render the same results as assuming no unobserved confounding.
#' @param out.formula Formula for the outcome regression model(s).
#' @param treat.formula Formula for the propensity score model (regression model for treatment assignment).
#' @param data data.frame containing the variables in the formulas.
#' @param Y Outcome vector. The arguments Y,T and X will be only used if out.formula, treat.formula or data matrix is not specified.
#' @param T Treatment assignment vector. The arguments Y,T and X will be only used will be used if out.formula, treat.formula or data matrix is not specified.
#' @param X Covariate matrix. The arguments Y,T and X will be only used if out.formula, treat.formula or data matrix is not specified.
#' @param rho0 Pre-specified value of \code{rho0}.
#' @param rho1 Pre-specified value of \code{rho1}.
#' @param param Specifies parameter of interest, it can be "ACE" (default value), "ACT", "Y0", "Y1", "Y0T1" or "Y1T0".
#' @param subset Specifies which subset of the covariates should be used in the regression models. It can be "noselection" (default value), "single", "double" or "refit". "Noselection" or "double" is suggested. See details.
#' @param sand Specifies which estimator of the standard errors should be used for OR, see details.
#' @param gridn Number of fixed points within the \code{rho0} and \code{rho1} intervals for which sigma0 and sigma1 should be estimated.
#' @param rho.plotrange  An interval larger than \code{rho0} and \code{rho1}. This specifies the range of value (on x axis) where the ui plot will be made using \code{\link{plot.ui}}.
#' @param alpha Default 0.05 corresponding to a confidence level of 95 for CI and UI.
#' @param sigma_correction A character variable which specifies if a corrected estimation of variance of the outcome model is used to find the confounding bias and if so whether we use the old or new correction method. It can take values 'non', 'old' or 'new'.
#' @param missing  This argument should not be changed manually. A logical variable which is FALSE by default and TRUE if this function is called withing the function \code{\link{ui.missing}}.
#'
#' @details In order to visualize the results, you can use \code{\link{plot.ui}}. Details about estimators can be found in Genbäck and de Luna (2018)
#'
#' The standard errors are calculated with the following estimators:
#'
#' DR - influence function-based estimator
#'
#' OR - if sand=TRUE sandwich estimator (default and recommended), if sand=FALSE large sample variance
#'
#' The set of covariates which are included in the outcome and propensity score models are choosed as below:
#'
#' noselection - All the covariates specified in formulas (and included in the data matrix) will be used in the respective models or the set of covariates in the X matrix will be used for all models.
#'
#' single - Lasso regressions will be fitted to the outcome models using covariates specified by formulas or X matrix and union of the selected sets will be used.
#'
#' double - Lasso regressions will be fitted to the outcome models and the propensity score model using covariates specified by formulas or X matrix and union of the selected sets will be used.
#'
#' refit- Lasso regressions will be fitted to the outcome models and the propensity score model using covariates specified by formulas or X matrix and union of the sets found by potential outcome models and the set found by the propensity score model will be used in refitting potential outcome and ps models with no regularizations respectively.
#'
#'
#' @return A list containing:
#' \item{call}{The matched call.}
#' \item{rho0}{The range of \code{rho0} from which the ui is calculated.}
#' \item{rho1}{If ACT==FALSE,range of \code{rho1} from which the ui is calculated.}
#' \item{out.model0}{Outcome regression model for non-treated.}
#' \item{out.model1}{Outcome regression model for treated.}
#' \item{treat.model}{Regression model for treatment mechanism (propensity score).}
#' \item{sigma0}{Estimate of sigma0 for different values of rho0.}
#' \item{sigma1}{Estimate of sigma1 for different values of rho1.}
#' \item{DR}{DR inference, confidence intervals for different pre-specified values of \code{rho0} and \code{rho1} for the OR estimator, uncertainty interval, coefficient estimates, confounding bias, identification interval, standard error etc.}
#' \item{OR}{OR inference, confidence intervals for different pre-specified values of \code{rho0} and \code{rho1} for the OR estimator, uncertainty interval, coefficient estimates, confounding bias, identification interval, standard error etc.}
#' \item{XYhat}{Union of the sets of covariates selected by fitting lasso regressions to the outcome models using using covariates specified by formulas or X matrix.}
#' \item{XThat}{Set of covariates selected by fitting lasso regression to the ps model using covariates specified by formulas or X matrix.}
#'
#' @references Genbäck, M., de Luna, X. (2018). Causal Inference Accounting for Unobserved Confounding after Outcome Regression and Doubly Robust Estimation. \emph{Biometrics}. DOI: 10.1111/biom.13001
#' @examples
#' library(MASS)
#' n <- 500
#' rho <- 0.3
#' x <- rnorm(n)
#' s0 <- 2
#' s1 <- 3
#' error <- mvrnorm(
#'   n, c(0, 0, 0),
#'   matrix(c(1, 0.6, 0.9, 0.6, 4, 0.54, 0.9, 0.54, 9), ncol = 3)
#' )
#' delta <- c(-0.3, 0.65)
#' tstar <- cbind(1, x) %*% delta + error[, 1]
#' t <- tstar > 0
#' y1 <- -2 + x - 3 * x^2 + error[, 3]
#' y0 <- 0.5 + 2 * x - x^2 + error[, 2]
#' y <- y0
#' y[t == 1] <- y1[t == 1]
#' data <- data.frame(y, t, x)
#'
#' mean(y1 - y0)
#' ui <- ui.causal(y ~ x + I(x^2), t ~ x + I(x^2),
#'   data = data, rho0 = c(0, 0.4), rho1 = c(0, 0.4),
#'   subset = "double", param = "ACE", sigma_correction = 'new'
#' )
#' ui <- ui.causal(
#'   Y = y, T = t, X = cbind(x, x^2),
#'   rho0 = c(0, 0.4), rho1 = c(0, 0.4),
#'   subset = "double", param = "ACE", sigma_correction = 'new'
#' )
#' ui
#' plot(ui)
#' ui$XYhat
#' ui$XThat
#'
#'
#' mean(y1[t == 1] - y0[t == 1])
#' ui <- ui.causal(
#'   Y = y, T = t, X = cbind(x, x^2),
#'   rho0 = c(0, 0.4), rho1 = c(0, 0.4),
#'   subset = "noselection", param = "ACT", sigma_correction = 'new'
#' )
#' ui
#' plot(ui)
#' @export

ui.causal <- function(out.formula = NULL, treat.formula = NULL,
                      data = NULL, Y = NULL, T = NULL, X = NULL,
                      rho0 = c(-0.1, 0.1), rho1 = c(-0.1, 0.1),
                      param = "ACE", subset = "noselection",
                      sand = TRUE, gridn = 21,
                      rho.plotrange = c(-0.5, 0.5), alpha = 0.05,
                      regularization_alpha = 1,
                      sigma_correction = 'new', missing = FALSE) {
  # comment this for now, one should be carefull when to use rho0 and when rho1, also it might be that we dont need ui in the plot (infer=F)
  # if (min(rho.plotrange) >= min(rho0)) {
  #   	warning('Lower bound of plotrange is >= lower bound of UI, it needs to be less than. You will not be able to plot this object.')
  #   rho.plotrange[1] <- (min(rho0) - 1) / 2
  # }
  # if (max(rho.plotrange) <= max(rho0)) {
  #   warning('Upper bound of plotrange is <= upper bound of UI, it needs to be greater than. You will not be able to plot this object.')
  #   rho.plotrange[2] <- (max(rho0) + 1) / 2
  # }

  y.data <- data.frame()
  t.data <- data.frame()
  Yname <- ""
  Tname <- ""
  XYnames_preselection <- c()
  XTnames_preselection <- c()
  out.formula_postselection <- NULL
  treat.formula_postselection <- NULL


  # out.formula=y~x+I(x^2);treat.formula=z~x+I(x^2);rho0=c(0,0.3);rho1=c(0,0.3); param="ACE"

  # Here we define data, and columnnames for XT XY Y and T, in case where the user chosen to 1)give formula and data or 2)X, Y and T seperately
  if (!is.null(out.formula) & !is.null(treat.formula) & !is.null(data)) {
    if (class(data) != "data.frame") {
      stop("Data must be a data frame")
    }

    names <- all.vars(out.formula)
    Yname <- names[1]
    Ymodelmatrix <- as.data.frame(model.matrix.lm(out.formula, na.action = NULL)[, -1, drop = FALSE])
    XYnames_preselection <- names(Ymodelmatrix)

    names <- all.vars(treat.formula)
    Tname <- names[1]
    Tmodelmatrix <- as.data.frame(model.matrix.lm(treat.formula, na.action = NULL)[, -1, drop = FALSE])
    XTnames_preselection <- names(Tmodelmatrix)

    covariatesdataframe <- cbind(Ymodelmatrix, Tmodelmatrix)
    data <- cbind(data[, Yname, drop = FALSE], data[, Tname, drop = FALSE], covariatesdataframe[, !duplicated(c(names(Ymodelmatrix), names(Tmodelmatrix))), drop = FALSE])
    names(data)[1:2] <- c(Yname, Tname)
  } else {
    if ((is.null(X) | is.null(Y) | is.null(T))) {
      stop("Variables cannot be NULL")
    }
    data <- data.frame(Y = Y, T = T, X)
    Yname <- "Y"
    Tname <- "T"
    XYnames_preselection <- XTnames_preselection <- names(data[-c(1, 2)])
  }
  # }

  if (subset != "noselection" & (length(XYnames_preselection) <= 1 | length(XTnames_preselection) <= 1)) {
    stop("Variable selection cannot be performed if the number of covariates is not bigger than one.")
  }

  # which indecies from XY and XT?
  XYhat <- XY1hat <- XY0hat <- XThat <- c()
  if (!subset == "noselection") {
    XY1hat <- vs.hdm(data[data[, Tname] == 1, Yname], data[data[, Tname] == 1, XYnames_preselection])
    if (!missing) {
        XY0hat <- vs.hdm(data[data[, Tname] == 0, Yname], data[data[, Tname] == 0, XYnames_preselection])
     }
    # XY1hat=vs.glmnet(data[data[,Tname]==1,Yname],data[data[,Tname]==1,XYnames_preselection],alpha=regularization_alpha)
    # XY0hat=vs.glmnet(data[data[,Tname]==0,Yname],data[data[,Tname]==0,XYnames_preselection],alpha=regularization_alpha)

    XYhat <- union(XY1hat, XY0hat)
  }
  if (!subset == "noselection" & !subset == "single") {
    XThat <- vs.glmnet(data[, Tname], data[, XTnames_preselection], alpha = regularization_alpha)
  }


  # if (subset != "noselection" & (length(XYhat) == 0 | length(XThat) == 0)) {
  #   stop("No covariates has been selected for at least one of the nuisance models.")
  # }

  outY_XY <- outT_XT <- ""
  if (subset == "noselection") {
    outY_XY <- XYnames_preselection
    outT_XT <- XTnames_preselection
  } else if (subset == "single") {
    if(length(XYhat) == 0) {
      stop("No covariates has been selected for the outcome models.")
    }
    outY_XY <- XYnames_preselection[XYhat]
    outT_XT <- XYnames_preselection[XYhat]
  } else if (subset == "double") {
   if( (length(XYhat) == 0 | length(XThat) == 0)) {
        stop("No covariates has been selected for at least one of the nuisance models.")
      }
    outY_XY <- union(XYnames_preselection[XYhat], XTnames_preselection[XThat])
    outT_XT <- union(XYnames_preselection[XYhat], XTnames_preselection[XThat])
  } else {
    if( (length(XYhat) == 0 | length(XThat) == 0)) {
      stop("No covariates has been selected for at least one of the nuisance models.")
    }
    outY_XY <- XYnames_preselection[XYhat]
    outT_XT <- XTnames_preselection[XThat]
  }

  FUN <- NULL
  if (param == "ACE") {
    FUN <- ui.ace
  } else if (param == "ACT") {
    FUN <- ui.act
  } else if (param == "Y0") {
    FUN <- ui.y0
  } else if (param == "Y1") {
    FUN <- ui.y1
  } else if (param == "Y0T1") {
    FUN <- ui.y0t1
  } else {
    FUN <- ui.y1t0
  }

  out.formula_postselection <- paste(Yname, paste(outY_XY, collapse = "+"), sep = " ~ ")
  treat.formula_postselection <- paste(Tname, paste(outT_XT, collapse = "+"), sep = " ~ ")


  # t.data<-get_all_vars(treat.formula_postselection,data=data)
  # y.data<-get_all_vars(out.formula_postselection,data=data)
  t.data <- data[, c(Tname, outT_XT)]
  y.data <- data[, c(Yname, outY_XY)]
  #  remove individuals with partial missing data in the covariates
  if (sum(complete.cases(y.data) == FALSE) > 0) {
    warning(paste("Partial missing values in covariates! ", sum(complete.cases(y.data) == FALSE), "individual(s) are removed from the outcome regression."))
    y.data <- y.data[complete.cases(y.data), , drop = FALSE]
  }
  if (sum(complete.cases(t.data) == FALSE) > 0) {
    warning(paste("Partial missing values in covariates. ", sum(complete.cases(t.data) == FALSE), "individual(s) are removed from the propensity score regression."))
    t.data <- t.data[complete.cases(t.data), , drop = FALSE]
  }



  treat.model <- glm(treat.formula_postselection, family = binomial(link = "probit"), data = t.data)
  # XThatdesign<-model.matrix(treat.model)
  gamma <- coef(summary(treat.model))[, 1]



  out.formula <- out.formula_postselection
  output <- FUN(out.formula, y.data, gamma, t.data,
    rho0 = rho0, rho1 = rho1,
    sand, gridn, rho.plotrange, alpha, sigma_correction
  )

  output$treat.model <- treat.model
  output$XYhat <- XYhat
  output$XThat <- XThat

  output$call <- match.call()
  class(output) <- "ui"


  return(output)
}
