#' Plot method for cuminc objects
#' 
#' @param x Object of class \code{"Cuminc"} to be printed or plotted
#' @param \dots Further arguments to plot or print method
#'   
#' @export 
plot.Cuminc <- function(x, ...)
{
  if (!inherits(x, "Cuminc"))
    stop("'x' must be a 'Cuminc' object")
  
  plot(attr(x, "survfit"), ...)
}
