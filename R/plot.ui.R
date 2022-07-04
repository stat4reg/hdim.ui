#' Plot ui objects
#'
#' Plot function for objects returned from \code{\link{ui.causal}}, \code{\link{ui.missing}}, \code{\link{ui.y0}}, \code{\link{ui.y1}}, \code{\link{ui.y0t1}}, \code{\link{ui.y1t0}}, \code{\link{ui.ace}} and \code{\link{ui.actT}}.
#' Plots confidence intervals for different values of rho and the uncertainty interval. Alternatively, a contour plot can be made in case of ACE.
#' @param x An object of class ui
#' @param DR If TRUE the doubly robust estimator is plotted, otherwise the outcome regression estimator is plotted.
#' @param main Main title, default is no title.
#' @param infer If TRUE confidence and uncertainty intervals are provided in chosen ranges of rho, default is TRUE.
#' @param lower A lower bound for the parameter of interest, if available.
#' @param upper An upper bound for the parameter of interest, if available.
#' @param contour If TRUE and ACE is the parameter of interest a contour plot will be made where rho0 and rho1 are on x and y axis, respectively. Default is FALSE.
#' @param ... Additional arguments, use is discouraged.
#'
#' @importFrom graphics arrows axis lines plot polygon text
#' @importFrom graphics contour
#' @export


plot.ui <- function(x, DR = TRUE, main = "", infer = TRUE,
                    lower = NULL, upper = NULL, contour = FALSE, xlab = NULL, ylab = NULL, ...) {


  # xlab=ylab=NULL
  if (DR == TRUE) {
    ex <- x$plot$DR
    # AA not used
    ui.obj <- x$DR
  } else {
    # AA DR to OR
    # ex<-x$plot$DR
    ex <- x$plot$OR
    # AA not used
    ui.obj <- x$OR
  }

  if (length(x$rho0) == 1 | length(x$rho1) == 1) {
    stop("Not able to plot ui object when rho takes only one value")
  }

  if (x$plot$plot == FALSE & contour == FALSE) {
    warning("Not able to plot ui object, with plot==FALSE. In ui.causal() rho0 and rho1 must be equal intervals otherwise plot is changed to FALSE.")
    gridrho0 <- x$plot$gridrho0[x$plot$nui0]
    gridrho1 <- x$plot$gridrho1[x$plot$nui1]
    contour(
      y = gridrho0, x = gridrho1, z = ui.obj$coef,
      ylab = expression(rho[0]), xlab = expression(rho[1]), main = main
    )
  } else {
    if (length(dim(ex$ci)) > 2 & contour == TRUE) {
      gridrho0 <- x$plot$gridrho0[x$plot$nui0]
      gridrho1 <- x$plot$gridrho1[x$plot$nui1]
      contour(
        y = gridrho0, x = gridrho1, z = ui.obj$coef,
        ylab = expression(rho[0]), xlab = expression(rho[1]), main = main
      )
    } else {
      if (contour == TRUE) {
        warning("Not able to plot contours for this parameter of interst.")
      }
      ## A if we never look at gridrho0 why not consider only gridrho in ui.causal
      # rho<-x$plot$gridrho0
      # define instead x$plot$gridrho based on which param is interested
      rho <- x$plot$gridrho

      ## A this never returns true since we fill the null rho values in the
      ## first lines of ui.causal
      # if(is.null(x$rho1)){
      if (length(dim(ex$ci)) == 2) {
        coef <- ex$coef
        ci <- cbind(ex$ci[, 1], ex$ci[, 2])
      } else {
        coef <- diag(ex$coef)
        ci <- cbind(diag(ex$ci[, , 1]), diag(ex$ci[, , 2]))
      }

      if (is.null(xlab)) {
        xlab <- expression(rho)
      }

      plot(rho, coef, type = "l", mgp = c(2, 1, 0), ylim = c(min(ci[, 1]), max(ci[, 2])), ylab = ylab, xlab = xlab, xaxs = "i", main = main)
      if (infer == TRUE) {
        polygon(c(rho, rev(rho)), c(ci[, 2], rev(ci[, 1])), col = "grey90", border = NA)

        lines(rho, coef)
        lines(rho, ci[, 1], lty = 2)
        lines(rho, ci[, 2], lty = 2)
        lines(c(-1, 1), c(0, 0))


        Stn <- x$plot$gridn
        Fin <- x$plot$gridn * 2 - 1
        n0 <- which(round(rho, 10) == 0)
        if (length(n0) > 1) {
          n0 <- which(abs(rho) == min(abs(rho)))
        }

        if (n0 == Stn | n0 == Fin) {
          nv <- Stn + round((Fin - Stn) / 4) * 0:4
          if (n0 == Stn) {
            nv <- nv[-1]
          } else {
            nv <- nv[-5]
          }
        } else {
          nv <- sort(Stn + round((Fin - Stn) / 3) * 0:3)
        }

        minR <- which(ci[Stn:Fin, 1] == min(ci[Stn:Fin, 1])) + (Stn - 1)
        maxR <- which(ci[Stn:Fin, 2] == max(ci[Stn:Fin, 2])) + (Stn - 1)

        lines(c(rho[4], rho[minR]), c(ci[minR, 1], ci[minR, 1]), lty = 2, col = "blue")
        lines(c(rho[4], rho[maxR]), c(ci[maxR, 2], ci[maxR, 2]), lty = 2, col = "blue")
        arrows(x0 = rho[Stn], y0 = min(ci[, 1]) - 10, y1 = ci[Stn, 1], code = 3, length = 0, lty = 2, col = "blue")
        arrows(x0 = rho[Fin], y0 = min(ci[, 1]) - 10, y1 = ci[Fin, 1], code = 3, length = 0, lty = 2, col = "blue")



        arrows(x0 = 0, y0 = ci[n0, 1], y1 = ci[n0, 2], code = 3, length = 0.1, col = "red")
        for (i in nv) {
          arrows(x0 = rho[i], y0 = ci[i, 1], y1 = ci[i, 2], code = 3, length = 0.1, lty = 3, lwd = 1.5)
        }

        arrows(x0 = rho[4], y0 = min(ci[Stn:Fin, 1]), y1 = max(ci[Stn:Fin, 2]), code = 3, length = 0.1, col = "blue")


        text(mean(rho[2:3]), mean(ci), "UI", col = "blue")
        axis(1, at = round(c(rho[Stn], rho[Fin]), 3), col.axis = "blue")
      }
      if (!is.null(lower)) {
        abline(lower, 0, col = "red")
      }
      if (!is.null(upper)) {
        abline(upper, 0, col = "blue")
      }
    }
  }
}
