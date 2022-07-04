#' Print function for object of class ui
#'
#' Print function for objects returned from \code{\link{ui.causal}}, \code{\link{ui.missing}}, \code{\link{ui.y0}}, \code{\link{ui.y1}}, \code{\link{ui.y0t1}}, \code{\link{ui.y1t0}}, \code{\link{ui.ace}} and \code{\link{ui.act}}.
#' @param x An object of class ui
#' @param digits number of digits to be printed.
#' @param digitsci number of digits to be printed in the confidence interval.
#' @param digitsui number of digits to be printed in the uncertainty interval.
#' @param ... Additional arguments, use is discouraged.
#'
#' @export


print.ui <- function(x, digits = 3, digitsci = digits, digitsui = digits, ...) {
  ci <- c(interv.p(x$OR$naivci, digitsci), interv.p(x$DR$naivci, digitsci))
  ui <- c(interv.p(x$OR$ui, digitsui), interv.p(x$DR$ui, digitsui))
  Tab <- data.frame(Est = round(c(x$OR$naiv, x$DR$naiv), digits), ci = ci, ui = ui)
  rownames(Tab) <- c("OR", "DR")
  # cat("\nCall:\n", deparse(x$call), "\n\n\n", sep = "")
  cat("Confidence intervals (CI) derived assuming unconfoundedness (rho=0)", "\n", sep = "")

  if (is.null(x$rho1)) {
    if (length(x$rho0) > 1) {
      cat("Uncertainty intervals (UI) derived assuming ", min(x$rho0), "<=rho0<=", max(x$rho0), "\n", "\n", sep = "")
    } else {
      cat("Uncertainty intervals (UI) derived assuming rho0=", x$rho0, "\n", "\n", sep = "")
    }
  } else {
    if (length(x$rho0) > 1 & length(x$rho1) > 1) {
      cat("Uncertainty intervals (UI) derived assuming ", min(x$rho0), "<=rho0<=", max(x$rho0), " and ", min(x$rho1), "<=rho1<=",
        max(x$rho1), "\n", "\n",
        sep = ""
      )
    } else {
      if (length(x$rho0) > 1) {
        cat("Uncertainty intervals (UI) derived assuming ", min(x$rho0), "<=rho0<=", max(x$rho0), " and rho1=", x$rho1, "\n", "\n", sep = "")
      } else {
        if (length(x$rho1) > 1) {
          cat("Uncertainty intervals (UI) derived assuming rho0=", x$rho0, " and ", min(x$rho1), "<=rho1<=",
            max(x$rho1), "\n", "\n",
            sep = ""
          )
        } else {
          cat("Uncertainty intervals (UI) derived assuming rho0=", x$rho0, " and rho1=", x$rho1, "\n", "\n", sep = "")
        }
      }
    }
  }
  print(Tab)
}
