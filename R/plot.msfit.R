#' Plot method for an msfit object
#' 
#' Plot method for an object of class \code{"msfit"}. It plots the estimated
#' cumulative transition intensities in the multi-state model.
#' 
#' @param x Object of class \code{"msfit"}, containing estimated cumulative transition
#' intensities for all transitions in a multi-state model
#' @param type One of \code{"single"} (default) or \code{"separate"}; in case
#' of \code{"single"}, all estimated cumulative hazards are drawn in a single
#' plot, in case of \code{"separate"}, separate plots are shown for the
#' estimated transition intensities
#' @param cols A vector specifying colors for the different transitions;
#' default is 1:K (K no of transitions), when type=\code{"single"}, and 1
#' (black), when type=\code{"separate"}
#' @param xlab A title for the x-axis; default is \code{"Time"}
#' @param ylab A title for the y-axis; default is \code{"Cumulative hazard"}
#' @param ylim The y limits of the plot(s); if ylim is specified for
#' type="separate", then all plots use the same ylim for y limits
#' @param lwd The line width, see \code{\link{par}}; default is 1
#' @param lty The line type, see \code{\link{par}}; default is 1
#' @param legend Character vector of length equal to the number of transitions,
#' to be used in a legend; if missing, these will be taken from the row- and
#' column-names of the transition matrix contained in \code{x$trans}. Also used
#' as titles of plots for type=\code{"separate"}
#' @param legend.pos The position of the legend, see \code{\link{legend}};
#' default is \code{"topleft"}
#' @param bty The box type of the legend, see \code{\link{legend}}
#' @param use.ggplot Default FALSE, set TRUE for ggplot version of plot
#' @param xlim Limits of x axis, relevant if use_ggplot = T
#' @param scale_type "fixed", "free", "free_x" or "free_y", see scales argument
#' of facet_wrap(). Only relevant for use_ggplot = T.
#' @param conf.int Confidence level (\%) from 0-1 for the cumulative hazard, 
#' default is 0.95. Only relevant for use_ggplot = T
#' @param conf.type Type of confidence interval - either "log" or "plain" . See
#' function details of \code{plot.probtrans} for details
#' @param \dots Further arguments to plot
#' 
#' @return No return value
#' 
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}

#' @seealso \code{\link{msfit}}
#' @keywords hplot
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # data in wide format, for transition 1 this is dataset E1 of
#' # Therneau & Grambsch (2000)
#' tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'         dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'         x1=c(1,1,1,0,0,0),x2=c(6:1))
#' # data in long format using msprep
#' tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
#' 		data=tg,keep=c("x1","x2"),trans=tmat)
#' # events
#' events(tglong)
#' table(tglong$status,tglong$to,tglong$from)
#' # expanded covariates
#' tglong <- expand.covs(tglong,c("x1","x2"))
#' # Cox model with different covariate
#' cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
#' 	data=tglong,method="breslow")
#' summary(cx)
#' # new data, to check whether results are the same for transition 1 as
#' # those in appendix E.1 of Therneau & Grambsch (2000)
#' newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
#' msf <- msfit(cx,newdata,trans=tmat)
#' # standard plot
#' plot(msf)
#' # specifying line width, color, and legend
#' plot(msf,lwd=2,col=c("darkgreen","darkblue","darkred"),legend=c("1->2","1->3","2->3"))
#' # separate plots
#' par(mfrow=c(2,2))
#' plot(msf,type="separate",lwd=2)
#' par(mfrow=c(1,1))
#' 
#' # ggplot version - see vignette for details
#' library(ggplot2)
#' plot(msf, use.ggplot = TRUE)
#' 
#' @export 
plot.msfit <- function(x, 
                       type = c("single", "separate"), 
                       cols,
                       xlab = "Time",
                       ylab = "Cumulative hazard", 
                       ylim,
                       lwd, 
                       lty,
                       legend, 
                       legend.pos = "right", 
                       bty = "n", 
                       use.ggplot = FALSE,
                       
                       # Possible ggplot args here
                       xlim,
                       scale_type = "fixed",
                       conf.int = 0.95, 
                       conf.type = "none",
                       ...) {
  
  # Prelim
  type <- match.arg(type)
  if (!inherits(x, "msfit")) stop("'x' must be a 'msfit' object")
  
  # Use ggplot 
  if (use.ggplot) {
    
    # Check for ggplot2
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package ggplot2 needed for this function to work. Please install it.", call. = FALSE)
    }
    
    p <- ggplot.msfit(
      x = x,
      type = type,
      cols = cols,
      xlab = xlab,
      ylab = ylab,
      ylim = ylim,
      lwd = lwd,
      lty = lty,
      legend = legend,
      legend.pos = legend.pos,
      xlim = xlim,
      scale_type = scale_type,
      conf.int = conf.int,
      conf.type = conf.type
    )
    return(p)
  }
  
  # Base R plot
  msf1 <- x$Haz
  K <- max(msf1$trans)
  msft <- unique(msf1$time) # the time points
  nt <- length(msft)
  msfp <- matrix(msf1[,2], nrow=nt) # the cumulative hazards in matrix (nt x K)
  if (missing(legend))
    legend <- to.trans2(x$trans)$transname
  if (type=="single") {
    if (missing(cols)) cols <- 1:K
    if (missing(ylim)) ylim <- c(0, max(msfp))
    if (missing(lwd)) lwd <- 1
    if (missing(lty)) lty <- rep(1, K)
    plot(msft, msfp[,1], type="s", ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1], lwd=lwd,
         lty=lty[1], ...)
    if (K > 1)
      for (k in 2:K)
        lines(msft, msfp[,k], type="s", col=cols[k], lwd=lwd, lty=lty[k], ...)
    if (missing(legend.pos))
      legend("topleft", legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
    else
      legend(legend.pos[1], legend.pos[2], legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
  }
  else if (type=="separate") {
    if (missing(cols)) cols <- rep(1,K)
    if (missing(lwd)) lwd <- 1
    if (missing(lty)) lty <- 1
    for (k in 1:K) {
      if (missing(ylim))
        plot(msft, msfp[,k], type="s", xlab=xlab, ylab=ylab, col=cols[k], lwd=lwd, ...)
      else
        plot(msft, msfp[, k], type="s", ylim=ylim, xlab=xlab, ylab=ylab, col=cols[k],
             lwd=lwd, ...)
      title(main=legend[k])
    }
  }
  return(invisible())
}
