#' Uncertainty intervals for partial correlation
#'
#' This function is used to calculate uncertainty intervals for partial correlation between two variables adjusting for a set of other variables.
#' The two variables of interest might have missing data according to the three missing data mechanisms (see reference). Data for the variables in the adjustment set should be completely observed.
#' @param out.formula Regression formula where one variable of interest is the outcome and another variable of interest is the first predictor. Other predictors are variables in the adjustment set.
#' @param mis.formula Formula for missingness mechanism. If NULL the same covariates as in the outcome regression will be used.
#' @param data data.frame containing the variables in the formula.
#' @param rho The min and the max of the sensitivity parameter.
#' @param rho2 The min and the max of the sensitivity parameter  \code{rho2} in the model for missing in the predictor of interest within missing data mechanism C. The same as \code{rho} by default.}
#' @param alpha Significance level. Default 0.05 corresponding to a confidence level of 95% for CI and UI.
#' @param gridn The number of distinct points  for \code{rho} at which confidence intervals should be constructed. Default is 101.
#' @param gridn2 The number of distinct points for \code{rho2} at which confidence intervals should be constructed for missing data mechanism C. Default is 11.
#' @details In order to visualize the results, you can use \code{\link{plot.uipcor}}. In the initial publication (see reference), the sensitivity parameter is called \gamma,  while the function denote the sensitivity paramter as \code{rho} for the correspondence with other functions in the package.
#' @return A list containing:
#' \item{call}{The matched call}
#' \item{out.formula}{Regression formula where one variable of interest is the outcome and another variable of interest is the first predictor. Other predictors are variables in the adjustment set. }
#' \item{out.model}{A result of a call to \code{lm} for the specified regression. }
#' \item{mis.formula}{Regression formula for missingness mechanism(s).}
#' \item{mis.model}{A result of a call to \code{glm} for probit model (s) for the missingness mechanism(s).}
#' \item{rho}{The range of values for the sensitivity parameter, \code{rho}, for which an uncertainty interval is constructed. Default  is (-0.3; 0.3).}
#' \item{gridrho}{The values of the sensitivity parameter \code{rho} (parameters \code{rho1} and \code{rho2} for missing data mechanisms C) for which confidence intervals are constructed.}
#' \item{pcor,rho0}{Estimated partial correlation assuming that the sensitivity parameter(s) is 0.}
#' \item{pcor}{Estimated partial correlation for different values of the sensitivity parameter(s) in \code{gridrho}.}
#' \item{ident.bound}{Bounds for the estimated identification region for partial correlation. An interval from the lowest to the largest estimated partial correlation. }
#' \item{pcor.se}{Standard error for different values of the sensitivity parameter(s) in \code{gridrho}.}
#' \item{ci.rho0}{Confidence interval for partial correlation assuming that the sensitivity parameter(s) is 0.}
#' \item{ci}{Confidence intervals for partial correlation for different values of the sensitivity parameter(s) in \code{gridrho}.}
#' \item{ui}{Uncertainty interval for partial correlation. An interval from the minimum lower bound to the maximum upper bound of estimated confidence intervals.}
#' @author Tetiana Gorbach
#' @references Gorbach, T., de Luna, X. (2018). Inference for partial correlation when data are missing not at random. \emph{Statistics & Probability Letters}, 141, 82-89.
#' @examples
#' library(MASS)
#' n <- 1000
#' rho <- 0.1
#' error <- mvrnorm(n, mu = c(0, 0, 0), Sigma = diag(c(1.16, 0.028^2 * (1 - rho^2), 1)))
#' X3 <- rnorm(n, mean = 67, sd = 7)
#' X4 <- rbinom(n, 1, prob = 0.3)
#' X2 <- 2.313 - 0.042 * X3 - 0.216 * X4 + error[, 1]
#' X1 <- 1.092 + 0.01 * X2 - 0.002 * X3 - 0.006 * X4 + 0.028 * rho * error[, 3] + error[, 2]
#' Z <- as.numeric(2.708 + 0.548 * X2 - 0.036 * X3 - 0.042 * X4 + error[, 3] > 0)
#' X1[Z == 0] <- NA
#' data <- data.frame(X1, X2, X3, X4)
#' ui <- ui.pcor(
#'   out.formula = X1 ~ X2 + X3 + X4,
#'   data = data,
#'   rho = c(0.1, 0.5),
#'   alpha = 0.05,
#'   gridn = 10
#' )
#' ui
#' @importFrom all.names get_all_vars lm glm stats binomial coef  model.matrix qnorm as.formula   Matrix rankMatrix
#' @export

ui.pcor <- function(out.formula, data, rho = c(-0.3, 0.3), rho2 = rho, alpha = 0.05, gridn = 101, gridn2 = 11) {
  variables.out.formula <- all.names(out.formula, functions = F)
  name.outcome <- variables.out.formula[1]
  name.predictor <- variables.out.formula[2]
  name.adjustment.set <- variables.out.formula[-c(1:2)]


  data.to.use <- get_all_vars(out.formula, data = data)
  #  Check parameters -------------------------------------------------------
  if (class(data) != "data.frame") {
    stop("Data must be a data frame.")
  }
  if (ncol(data.to.use) == 2) {
    stop("No variables in the adjustment set according to the outcome formula.")
  }
  if (any(abs(rho) >= 1)) {
    stop("The limits for rho should be between -1 and 1.")
  }
  if (length(rho) > 2 | length(rho2) > 2) {
    stop("More than two numbers are used to specify the min and max of the sensitivity parameters.")
  }

  data.adjustment.set <- data.to.use[, name.adjustment.set]


  # Define Missing Data Mechanism -------------------------------------------
  # stop if some values for the variables in the adjustment set are missing
  if (any(is.na(data.adjustment.set))) {
    if (ncol(data.to.use) == 3) {
      variables.with.missing <- names(data.to.use)[-c(1:2)]
    } else {
      variables.with.missing <- names(data.adjustment.set[, colSums(is.na(data.adjustment.set)) > 0])
    }
    stop(paste0("The method is developed for the case of no missing in the variables in the adjustment set. Variables ",
      paste0(variables.with.missing, collapse = ", "),
      " have missing values. Suggestion: use data[complete.cases(data[, -c(1:2)]), ] as a data frame.",
      collapse = ""
    ))
  }
  # stop if not missing in the outcome, else MDM A if missing only in the outcome,
  # MDM B if the same missing in the outcome and the first predictor,
  # MDM C if different missing in the outocme and the first predictor
  if (all(!is.na(data.to.use[, name.outcome]))) {
    stop("No missing in the outcome")
  } else {
    if (all(!is.na(data.to.use[, 2]))) {
      print("There are missing values in the outcome, but no missingness in the predictor of interest. Calculating estimates according to missing data mechanism A.")
      mdm <- "A"
    } else {
      z1 <- is.na(data.to.use[, 1])
      z2 <- is.na(data.to.use[, 2])
      if (all(z1 == z2)) {
        print("The same missingness in both variables of interest, calculating estimates according to missing data mechanism B.")
        mdm <- "B"
      } else {
        if (any(z1 == z2)) {
          print("Different missingness in the two variables of interest, calculating estimates according to missing data mechanism C.")
          mdm <- "C"
        } else {
          stop("No observations with both variables of interest observed.")
        }
      }
    }
  }

  # main --------------------------------------------------------------------
  output <- list()


  if (mdm == "A" || mdm == "B") {
    # define gridrho
    if (length(rho) > 1) {
      gridrho <- rho[1] + ((rho[2] - rho[1]) / (gridn - 1)) * 0:(gridn - 1)
    } else {
      gridrho <- rho
      gridn <- 1
    }
    gridrho <- c(0, gridrho)
    # calculate estimates
    out.model <- lm(out.formula, data = data.to.use[!is.na(data.to.use[, name.outcome]), ])

    z <- as.numeric(!is.na(data.to.use[, name.outcome]))
    data$z <- z
    if (mdm == "A") {
      mis.formula <- update.formula(out.formula, as.formula(paste("z", "~ .")))
    } else {
      mis.formula <- as.formula(paste0("z~", paste0(name.adjustment.set, collapse = "+")))
    }
    mis.model <- glm(mis.formula, family = binomial(link = "probit"), data = data)
    output$mis.formula <- mis.formula
    output$mis.model <- mis.model


    X <- model.matrix(out.model)
    Xz <- model.matrix(mis.model)[z == 1, ]

    # estimating delta and u
    delta <- coef(summary(mis.model))[, 1]
    u <- (Xz %*% delta) # similar to other functions in the package

    # Estimating BetaOLS
    BetaOLS <- coef(summary(out.model))[2, 1]
    sigmaOLS <- summary(out.model)$sigma

    # Estimation of the term with lambda1 for bound
    OLSlambda <- solve(t(X) %*% X) %*% t(X) %*% lambda1(u)
    OLSlambda <- OLSlambda[2]
    sigma1.2p <- sigmaOLS /
      sqrt(1 + gridrho^2 * c(
        (-t(u) %*% lambda1(u) - t(lambda1(u)) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% lambda1(u)) /
          (dim(X)[1] - dim(X)[2])
      )) # sigmaOLScor1(X, sigmaOLS, n = dim(X)[1], p = dim(X)[2] , u, gridrho)
    if (length(rho) > 1) {
      Mt <- matrix(gridrho * sigma1.2p, ncol = 1) %*% matrix(OLSlambda, nrow = 1)
      Beta2 <- apply(Mt, 1, function(x) BetaOLS - x)
    } else {
      Beta2 <- BetaOLS - gridrho * as.vector(sigma1.2p) * OLSlambda
    }
    sigma2.3p <- summary(lm(as.formula(paste0(name.predictor, "~", paste0(name.adjustment.set, collapse = "+"))),
      data = data
    ))$sigma
    pcor <- Beta2 / sqrt(Beta2^2 + sigma1.2p^2 / sigma2.3p^2)


    pcor.se <- sqrt(sigma1.2p^2 *
      (1 - gridrho^2 * c(t(u) %*% lambda1(u) - t(lambda1(u)) %*% lambda1(u)) / dim(X)[1]) *
      solve(t(X) %*% X)[2, 2] /
      (Beta2^2 + sigma1.2p^2 / sigma2.3p^2))
    gridrho <- gridrho[-1]
    output$ci <- array(dim = c(2, length(gridrho)), dimnames = list(c("lower", "upper"), rho = as.character(round(gridrho, 4))))
  }

  if (mdm == "C") {
    gridrho <- matrix(t(
      expand.grid(
        rho[1] + ((rho[2] - rho[1]) / (gridn - 1)) * 0:(gridn - 1),
        rho2[1] + ((rho2[2] - rho2[1]) / (gridn2 - 1)) * 0:(gridn2 - 1)
      )
    ),
    nrow = 2,
    dimnames = list(c("rho1", "rho2"), NULL)
    )
    gridrho <- cbind(c(0, 0), gridrho)
    out.model <- lm(out.formula, data = data.to.use[!is.na(data.to.use[, c(name.outcome, name.predictor)]), ])

    z1 <- as.numeric(!is.na(data.to.use[, name.outcome])) # define the variable for missing outcome
    data$z1 <- z1
    z2 <- as.numeric(!is.na(data.to.use[, name.predictor])) # define the variable for missing predictor
    data$z2 <- z2

    mis.formula1 <- as.formula(paste0("z1~", paste0(name.adjustment.set, collapse = "+")))
    mis.formula2 <- as.formula(paste0("z2~", paste0(name.adjustment.set, collapse = "+")))
    mis.model1 <- glm(mis.formula1, family = binomial(link = "probit"), data = data)
    mis.model2 <- glm(mis.formula2, family = binomial(link = "probit"), data = data)

    output$mis.formula <- list(mis.formula1, mis.formula2)
    output$mis.model <- list(mis.model1, mis.model2)

    X <- model.matrix(out.model)
    Xz1 <- model.matrix(mis.model1)[z1 == 1, ]
    Xz2 <- model.matrix(mis.model2)[z2 == 1, ]
    X.complete.cases <- model.matrix(mis.model1)[z1 == 1 & z2 == 1, ]

    # estimating delta and u
    delta1 <- coef(summary(mis.model1))[, 1]
    u <- (X.complete.cases %*% delta1)

    delta2 <- coef(summary(mis.model2))[, 1]
    w <- (Xz2 %*% delta2)

    # Estimating BetaOLS
    BetaOLS <- coef(summary(out.model))[2, 1]
    sigmaOLS <- summary(out.model)$sigma

    # Estimation of the term with lambda1 for bound
    OLSlambda <- solve(t(X) %*% X) %*% t(X) %*% lambda1(u)
    OLSlambda <- OLSlambda[2]

    sigma1.2p <- sigmaOLS /
      sqrt(1 + gridrho[1, ]^2 * c(
        (-t(u) %*% lambda1(u) - t(lambda1(u)) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% lambda1(u)) /
          (dim(X)[1] - dim(X)[2])
      )) # sigmaOLScor1(X, sigmaOLS, n = dim(X)[1], p = dim(X)[2] , u, gridrho[1, ])

    if (dim(gridrho)[2] > 1) {
      Mt <- matrix(gridrho[1, ] * sigma1.2p, ncol = 1) %*% matrix(OLSlambda, nrow = 1)
      Beta2 <- apply(Mt, 1, function(x) BetaOLS - x)
    } else {
      Beta2 <- BetaOLS - gridrho[1, 1] * as.vector(sigma1.2p) * OLSlambda
    }

    sigma2.3p <- summary(lm(as.formula(paste0(name.predictor, "~", paste0(name.adjustment.set, collapse = "+"))),
      data = data
    ))$sigma /
      sqrt(1 + gridrho[2, ]^2 *
        c((-t(w) %*% lambda1(w) - t(lambda1(w)) %*% Xz2 %*% solve(t(Xz2) %*% Xz2) %*% t(Xz2) %*% lambda1(w)) /
          (dim(Xz2)[1] - dim(Xz2)[2])))
    pcor <- Beta2 / sqrt(Beta2^2 + sigma1.2p^2 / sigma2.3p^2)
    pcor.se <- sqrt(sigma1.2p^2 *
      (1 - gridrho[1, ]^2 * c(t(u) %*% lambda1(u) - t(lambda1(u)) %*% lambda1(u)) / dim(X)[1]) *
      solve(t(X) %*% X)[2, 2] /
      (Beta2^2 + sigma1.2p^2 / sigma2.3p^2))
    gridrho <- gridrho[, -1]
    output$ci <- array(dim = c(2, dim(gridrho)[2]), dimnames = list(c("lower", "upper"), NULL))
  }



  output$call <- match.call()
  output$out.formula <- out.formula
  output$out.model <- out.model

  output$rho <- rho
  output$rho2 <- rho2
  output$gridrho <- gridrho

  output$pcor.rho0 <- pcor[1]
  output$pcor.se.rho0 <- pcor.se[1]

  output$pcor <- pcor[-1]
  output$ident.bound <- cbind(min(output$pcor), max(output$pcor))

  output$pcor.se <- pcor.se[-1]
  zalpha <- qnorm(1 - alpha / 2)
  output$ci.rho0 <- c(output$pcor.rho0 - zalpha * output$pcor.se.rho0, output$pcor.rho0 + zalpha * output$pcor.se.rho0)
  output$ci[1, ] <- output$pcor - zalpha * output$pcor.se
  output$ci[2, ] <- output$pcor + zalpha * output$pcor.se
  output$ui <- cbind(min(output$ci[1, ]), max(output$ci[2, ]))

  colnames(output$ident.bound) <- c("lower", "upper")
  colnames(output$ui) <- c("lower", "upper")
  output <- output[c("call", "out.formula", "out.model", "mis.formula", "mis.model", "rho", "gridrho", "pcor.rho0", "pcor", "ident.bound", "pcor.se", "ci.rho0", "ci", "ui")]
  class(output) <- "uipcor"
  output
}
