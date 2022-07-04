#' Prints objects of class uipcor
#'
#' @param x an objects returned from \code{\link{ui.pcor}}
#' @param digits number of digits to be printed.
#' @param digitsci number of digits to be printed in the confidence interval.
#' @param digitsui number of digits to be printed in the uncertainty interval.
#' @param ... Additional arguments, use is discouraged.
#'
#' @export

print.uipcor <- function(x, digits = 3, digitsci = digits, digitsui = digits, ...) {
  ui <- apply(x$ui, 1, function(x) {
    interv.p(x, digitsui)
  })
  ci <- interv.p(x$ci.rho0, digitsci)

  cat("\nCall:\n", deparse(x$call), "\n\n\n", sep = "")
  cat("The estimate of the partial correlation and the printed below confidence interval (CI) are derived assuming ignorable dropout (rho=0)", "\n", sep = "")
  cat("Uncertainty intervals (UI) derived assuming ", x$rho[1], "<=rho<=", x$rho[2], "\n", "\n", sep = "")
  Tab <- data.frame(Est = round(x$pcor.rho0, digits), ci = ci, ui = ui)
  print(Tab)
}
