#' Plot of UI and CI for partial correlation
#'
#' Plot function for objects returned from \code{\link{ui.pcro}}. Plots confidence intervals,
#' coefficients and significans assuming ignorability and the uncertainty interval under non-ignorability.
#' @param x An object of class uiols
#' @param col Vector containing the color of confidence intervals (default black) and uncertainty intervals (default red).
#' @param ... Additional arguments, use is discouraged.
#'
#'
#' @importFrom graphics arrows axis lines par plot points text
#' @importFrom stats pnorm
#' @export

plot.uipcor <- function(x, plot.all = TRUE, which = NA, intercept = FALSE, ylab = NULL, col = c("black", "red"), ...) {
  ci <- x$ci.rho0
  ui <- x$ui

  cexM <- 1
  plot(ui, c(1, 1),
    type = "l",
    xlim = c(min(ci, ui), max(ci, ui)),
    ylim = c(1 - 0.2, 1 + 0.1),
    col = col[1],
    lwd = 2 * cexM,
    cex.main = cexM,
    cex.axis = cexM,
    xlab = "",
    ylab = "",
    yaxt = "n",
    main = ""
  )
  lines(ci, c(0.9, 0.9),
    col = col[2],
    lwd = 2 * cexM
  )
  points(x = x$pcor.rho0, y = 0.9, pch = 19, cex = 0.6 * cexM, col = col[2])
  axis(2, at = c(1, 0.9), labels = c("UI", "CI"), las = 1, cex.axis = cexM)
}
