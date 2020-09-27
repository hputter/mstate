#' Plot method for cuminc objects
#' 
#' @param x Object of class \code{"Cuminc"} to be printed or plotted
#' @inheritParams plot.probtrans
#' @param legend Character vector corresponding to number of absorbing states.
#' In case of a grouped Cuminc() object, with facet = F the length of the vector
#' is number absorbing states * group levels. Only relevant when use.ggplot = T
#' @param facet In case of group used for Cuminc, facet by it - 
#' only relevant when use.ggplot = T
#' @param cols Vector (numeric or character) specifying colours of the lines
#' @param \dots Further arguments to plot or print method
#'   
#' @export 
plot.Cuminc <- function(x,
                        use.ggplot = F,
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
                        facet = F,
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
