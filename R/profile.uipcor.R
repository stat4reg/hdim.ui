#' Plot of confidence intervals for partial correlation by values of the sensitivity parameter (s).
#'
#' Plots confidence intervals by values of sensitivity parameter (s) for objects returned from \code{\link{ui.pcor}}.
#' @param x An object of class uipcor
#' @param ... Additional arguments, for instance margins.
#' @details Plots estimated partial correlation, lower and upper bounds of confidence intervals for different values of rho (based on on values of the sensitivity parameter(s) in \code{x$gridrho}).
#' Plots a 2D plot when missing data only in the outcome or the same missing in the outcome and the predictor of interest (missing data mechanisms A and B).
#' Plots a 3D plot when different missing data in the outcome and the predictor of interest (missing data mechanism C).
#'
#' @importFrom dplyr graphics lines par plot polygon plotly
#' @export


profile.uipcor <- function(x, ...) {
  if (!is.matrix(x$gridrho)) {
    plot(x$gridrho, x$pcor,
      type = "l",
      ylim = c(min(x$ui), max(x$ui)),
      xlab = expression(rho), ylab = "Estimated partial correlation"
    )
    polygon(c(x$gridrho, rev(x$gridrho)), c(x$ci[2, ], rev(x$ci[1, ])), col = "grey90", border = NA)
    lines(x$gridrho, x$pcor)
    lines(x$gridrho, x$ci[1, ], lty = 2)
    lines(x$gridrho, x$ci[2, ], lty = 2)
  } else {
    ci_lower <- matrix(x$ci[1, ], ncol = length(unique(x$gridrho[1, ])), byrow = T)
    ci_upper <- matrix(x$ci[2, ], ncol = length(unique(x$gridrho[1, ])), byrow = T)

    plot_ly(
      colors = c("red", "blue"),
      x = unique(x$gridrho[1, ]),
      y = unique(x$gridrho[2, ]),
      z = matrix(x$pcor, ncol = length(unique(x$gridrho[1, ])), byrow = T)
    ) %>%
      add_surface(showscale = F) %>%
      add_surface(z = ~ci_lower, opacity = 0.98, colorscale = list(c(0, 1), c("rgb(128,128,128)", "rgb(128,128,120)")), showscale = F) %>%
      add_surface(z = ~ci_upper, opacity = 0.98, colorscale = list(c(0, 1), c("rgb(128,128,128)", "rgb(128,128,120)")), showscale = F) %>%
      layout(scene = list(
        xaxis = list(title = "rho_1"),
        yaxis = list(title = "rho_2"),
        zaxis = list(title = "Estimated partial correlation")
      ))
  }
}
