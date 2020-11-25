#' Plot method for Cuminc objects
#' 
#' Plot the estimates of the non-parametric Aalen-Johansen estimate of the 
#' cumulative incidence functions (competing risks data). Note this is a method
#' for \code{mstate::Cuminc} and not \code{cmprsk::cuminc}. Both return the same 
#' estimates, though the former does so in a dataframe, and the latter in the list.
#' 
#' Grouped cumulative incidences can be plotted either in the same plot or in facets,
#' see the \code{facet} argument.
#' 
#' @param x Object of class \code{"Cuminc"} to be printed or plotted
#' @inheritParams plot.probtrans
#' @param legend Character vector corresponding to number of absorbing states.
#' In case of a grouped \code{"Cuminc"} object, with facet = FALSE the 
#' length of the vector is number absorbing states * group levels. 
#' Only relevant when use.ggplot = TRUE
#' @param facet Logical, in case of group used for \code{"Cuminc"}, facet by it - 
#' only relevant when use.ggplot = TRUE
#' @param cols Vector (numeric or character) specifying colours of the lines
#' @param \dots Further arguments to plot or print method
#' 
#' @return A ggplot object if use.ggplot = T used, otherwise NULL.
#' 
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#' 
#' @examples 
#' library(ggplot2)
#' 
#' data("aidssi")
#' head(aidssi)
#' si <- aidssi
#' 
#' # No grouping
#' cum_incid <- Cuminc(
#' time = "time",
#' status = "status",
#' data = si
#' )
#' 
#' plot(
#' x = cum_incid,
#' use.ggplot = TRUE,
#' conf.type = "none",
#' lty = 1:2,
#' conf.int = 0.95
#' )
#' 
#' # With grouping
#' cum_incid_grp <- Cuminc(
#' time = "time",
#' status = "status",
#' group = "ccr5",
#' data = si
#' )
#' 
#' plot(
#'  x = cum_incid_grp,
#'  use.ggplot = TRUE,
#'  conf.type = "none",
#'  lty = 1:4, 
#'  facet = TRUE 
#' )
#'   
#' @export 
plot.Cuminc <- function(x,
                        use.ggplot = FALSE,
                        xlab = "Time",
                        ylab = "Probability",
                        xlim,
                        ylim,
                        lty,
                        legend,
                        cols, 
                        conf.type = c("log", "plain", "none"),
                        conf.int = 0.95,
                        legend.pos = "right",
                        facet = FALSE,
                        ...)
{
  if (!inherits(x, "Cuminc"))
    stop("'x' must be a 'Cuminc' object")
  
  # Ggplot version
  if (use.ggplot) {
    conf.type <- match.arg(conf.type)
    
    p <- ggplot.Cuminc(
      x = x,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = ylim,
      lty = lty,
      legend = legend,
      cols = cols, 
      conf.type = conf.type,
      conf.int = conf.int,
      legend.pos = legend.pos,
      facet = facet
    )
    
    return(p)
  } else {
    
    # Base r
    plot(attr(x, "survfit"), ...)
  }
}
