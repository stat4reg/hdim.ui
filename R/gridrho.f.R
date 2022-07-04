#' Divides an interval into a grid
#'
#' This is a support function for \code{\link{ui.y0}}, \code{\link{ui.y1}}, \code{\link{ui.y0t1}} and \code{\link{ui.y1t0}} that divides the \code{rho} interval into a grid
#' @param rho interval that should be divided
#' @param gridn number of grid points
#' @param rho.plotrange a larger interval of grids to be used in a plot
#' @param plot whether or not the larger interval of grids should be created
#'
#' @export


gridrho.f <- function(rho, gridn, rho.plotrange, plot) {
  if (length(rho) > 1) {
    gridrho <- rho[1] + ((rho[2] - rho[1]) / (gridn - 1)) * (0:(gridn - 1))
    if (plot == TRUE) {
      gridrho <- c(rho.plotrange[1] + ((rho[1] - rho.plotrange[1]) / (gridn - 1)) * (0:(gridn - 2)), gridrho)
      gridrho <- c(gridrho, rho[2] + ((rho.plotrange[2] - rho[2]) / (gridn - 1)) * (1:(gridn - 1)))
    }
    # moved here from below, so that ACE with rho.length=1 does not return error
    if (sum(gridrho == 0) == 0) {
      gridrho <- sort(c(gridrho, 0))
    }
  } else {
    gridrho <- rho
    # gridrho<-paste("rho0=",rho0," and ","rho1=",rho1,sep="")
  }
  # if(sum(gridrho==0)==0){
  # 	gridrho<-sort(c(gridrho,0))
  # }
  # vector with logical values of which values of rho should be included in the UI
  nui <- gridrho >= min(rho) & gridrho <= max(rho)
  return(list(gridrho, nui))
}
