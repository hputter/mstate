fillplot <- function(x,y1,y2,col,lwd)
  # y2>y1, x (ascending order), y1, y2 same length, added to existing plot, intended for type="s"
{
  nx <- length(x)
  # add mini-bit of space, this is to incorporate the possibility of a jump at the end 
  x <- c(x, x[nx]+0.1*diff(range(x))) 
  xx <- c(rep(x,c(1,rep(2,nx-1),1)),rep(rev(x),c(1,rep(2,nx-1),1)))
  yy <- c(rep(y1,rep(2,nx)),rep(rev(y2),rep(2,nx)))
  polygon(xx,yy,col=col,lwd=lwd)
}





#' Plot method for a probtrans object
#' 
#' Plot method for an object of class 'probtrans'. It plots the transition
#' probabilities as estimated by \code{\link{probtrans}}.
#' 
#' Regarding confidence intervals: let \eqn{p} denote a predicted probability,
#' \eqn{\sigma} its estimated standard error,
#' and \eqn{z_{\alpha/2}} denote the critical value of the standard normal 
#' distribution at confidence level \eqn{1 - \alpha}.
#' 
#' The confidence interval of type "plain" is then
#' \deqn{p \pm z_{\alpha/2} * \sigma}
#' 
#' The confidence interval of type "log", based on the Delta method, is then
#' \deqn{\exp(\log(p) \pm z_{\alpha/2} * \sigma / p)}
#' 
#' @param x Object of class 'probtrans', containing estimated transition
#' probabilities
#' @param from The starting state from which the probabilities are used to plot
#' @param type One of \code{"stacked"} (default), \code{"filled"},
#' \code{"single"} or \code{"separate"}; in case of \code{"stacked"}, the
#' transition probabilities are stacked and the distance between two adjacent
#' curves indicates the probability, this is also true for \code{"filled"}, but
#' the space between adjacent curves are filled, in case of \code{"single"},
#' the probabilities are shown as different curves in a single plot, in case of
#' \code{"separate"}, separate plots are shown for the estimated transition
#' probabilities
#' @param ord A vector of length equal to the number of states, specifying the
#' order of plotting in case type=\code{"stacked"} or \code{"filled"}
#' @param cols A vector specifying colors for the different transitions;
#' default is a palette from green to red, when type=\code{"filled"} (reordered
#' according to \code{ord}, and 1 (black), otherwise
#' @param xlab A title for the x-axis; default is \code{"Time"}
#' @param ylab A title for the y-axis; default is \code{"Probability"}
#' @param xlim The x limits of the plot(s), default is range of time
#' @param ylim The y limits of the plot(s); if ylim is specified for
#' type="separate", then all plots use the same ylim for y limits
#' @param lwd The line width, see \code{\link{par}}; default is 1
#' @param lty The line type, see \code{\link{par}}; default is 1
#' @param cex Character size, used in text; only used when
#' type=\code{"stacked"} or \code{"filled"}
#' @param legend Character vector of length equal to the number of transitions,
#' to be used in a legend; if missing, numbers will be used; this and the
#' legend arguments following are ignored when type="separate"
#' @param legend.pos The position of the legend, see \code{\link{legend}};
#' default is \code{"topleft"}
#' @param bty The box type of the legend, see \code{\link{legend}}
#' @param xaxs See \code{\link{par}}, default is "i", for type=\code{"stacked"}
#' @param yaxs See \code{\link{par}}, default is "i", for type=\code{"stacked"}
#' @param use.ggplot Default FALSE, set TRUE for ggplot version of plot
#' @param conf.int Confidence level (\%) from 0-1 for probabilities, 
#' default is 0.95 (95\% CI). Setting to 0 removes the CIs.
#' @param conf.type Type of confidence interval - either "log" or "plain" . See
#' function details for details. 
#' @param label Only relevant for type = "filled" or "stacked", set to 
#' "annotate" to have state labels on plot, or leave unspecified.
#' @param \dots Further arguments to plot
#' 
#' @return No return value
#' 
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#' 
#' @seealso \code{\link{probtrans}}
#' @keywords hplot
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # data in wide format, for transition 1 this is dataset E1 of
#' # Therneau and Grambsch (2000)
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
#' # those in appendix E.1 of Therneau and Grambsch (2000)
#' newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
#' msf <- msfit(cx,newdata,trans=tmat)
#' # probtrans
#' pt <- probtrans(msf,predt=0)
#' # default plot
#' plot(pt,ord=c(2,3,1),lwd=2,cex=0.75)
#' # filled plot
#' plot(pt,type="filled",ord=c(2,3,1),lwd=2,cex=0.75)
#' # single plot
#' plot(pt,type="single",lwd=2,col=rep(1,3),lty=1:3,legend.pos=c(8,1))
#' # separate plots
#' par(mfrow=c(2,2))
#' plot(pt,type="sep",lwd=2)
#' par(mfrow=c(1,1))
#' 
#' # ggplot version - see vignette for details
#' library(ggplot2)
#' plot(pt, ord=c(2,3,1), use.ggplot = TRUE)
#' 
#' @export 
plot.probtrans <- function(x, 
                           from = 1, 
                           type = c("filled", "single", "separate", "stacked"),
                           ord, # ord for "stacked" and "filled "only, cex for text only
                           cols,
                           xlab = "Time", 
                           ylab = "Probability", 
                           xlim, 
                           ylim, 
                           lwd, 
                           lty, 
                           cex, 
                           legend, 
                           legend.pos = "right",
                           bty = "n", 
                           xaxs = "i", 
                           yaxs = "i", 
                           use.ggplot = FALSE,
                           
                           # Ggplot args here
                           conf.int = 0.95,
                           conf.type = c("log", "plain", "none"),
                           label, # "annotate" if want like base R
                           ...) {
  
  # Prelims
  type <- match.arg(type)
  if (!inherits(x, "probtrans"))
    stop("'x' must be a 'probtrans' object")
  
  # Ggplot version
  if (use.ggplot) {
    
    # Check for ggplot2
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package ggplot2 needed for this function to work. Please install it.", call. = FALSE)
    }
    
    conf.type <- match.arg(conf.type)
    
    p <- ggplot.probtrans(
      x = x,
      from = from,
      type = type, 
      ord = ord,
      cols = cols,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = ylim,
      lwd = lwd, # linewidth
      lty = lty,
      cex = cex, # size of label
      legend = legend,
      legend.pos = legend.pos,
      conf.int = conf.int, # 0 if no confidence intervals
      conf.type = conf.type,
      label = label
    )
    
    return(p)
  }
  
  trans <- x$trans
  S <- dim(trans)[1]
  if ((from<1) | (from>S)) stop("'from' incorrect")
  pt1 <- x[[from]]
  ptt <- pt1$time # the time points
  nt <- length(ptt)
  ptp <- pt1[,2:(S+1)] # those are the actual transition probabilities
  if (missing(ord)) ord <- 1:S
  if (missing(lwd)) lwd <- 1
  lwd <- rep(lwd, S)
  lwd <- lwd[1:S]
  if (missing(lty)) lty <- 1
  lty <- rep(lty, S)
  lty <- lty[1:S]
  if (missing(cex)) cex <- 1
  cex <- rep(cex, S)
  cex <- cex[1:S]
  if (missing(cols)) {
    if (type=="filled" | type=="single") cols <- rev(RColorBrewer::brewer.pal(S, "RdYlGn"))
    else cols <- 1
  }
  cols <- rep(cols, S)
  cols <- cols[1:S]
  if (missing(legend)) {
    legend <- dimnames(trans)[[2]]
    if (is.null(legend)) legend <- as.character(1:S)
  }
  else if (length(legend) != S) stop("legend has incorrect length")
  if (type=="single") {
    if (missing(xlim)) xlim <- range(ptt)
    if (missing(ylim)) ylim <- c(0,max(ptp))
    plot(ptt, ptp[,1], type="s", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1],
         lwd=lwd[1], lty=lty[1], ...)
    for (s in 2:S) lines(ptt, ptp[,s], type="s", col=cols[s], lwd=lwd[s], lty=lty[s], ...)
    if (missing(legend.pos))
      legend("topright", legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
    else
      legend(legend.pos[1], legend.pos[2], legend=legend, col=cols, lwd=lwd, lty=lty, bty=bty)
  }
  else if (type=="stacked") {
    if (missing(xlim)) xlim <- range(ptt)
    if (missing(ylim)) ylim <- c(0,1)
    # finally change order according to ord
    lwd <- lwd[ord]; lty <- lty[ord]; cex <- cex[ord]; cols <- cols[ord]
    y0 <- 0
    ptpsum <- ptp[,ord[1]]
    eps <- (xlim[2]-xlim[1])/50
    nt <- sum(ptt <= xlim[2]-eps)
    dy <- ptp[nt,ord[1]]
    y <- y0 + dy/2
    y1 <- y0 + dy
    plot(ptt, ptpsum, type="s", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1],
         lwd=lwd[1], ...)
    text(xlim[2]-eps, y, legend[ord[1]], adj=1, cex=cex)
    for (s in 2:S) {
      ptpsum <- ptpsum + ptp[,ord[s]]
      lines(ptt, ptpsum, type="s", col=cols[s], lwd=lwd[s], ...)
      y0 <- y1
      dy <- ptp[nt,ord[s]]
      y <- y0 + dy/2
      y1 <- y0 + dy
      text(xlim[2]-eps, y, legend[ord[s]], adj=1, cex=cex)
    }
  }
  else if (type=="filled") {
    par(xaxs=xaxs, yaxs=yaxs)
    if (missing(xlim)) xlim <- range(ptt)
    if (missing(ylim)) ylim <- c(0,1)
    # finally change order according to ord
    lwd <- lwd[ord]; lty <- lty[ord]; cex <- cex[ord]; cols <- cols[ord]
    y0 <- 0
    eps <- (xlim[2]-xlim[1])/50
    nt <- sum(ptt <= xlim[2])
    ptt <- ptt[1:nt]
    ptplow <- rep(0,nt)
    ptpup <- ptp[1:nt, ord[1]]
    dy <- ptp[nt, ord[1]]
    y <- y0 + 0.5*dy
    y1 <- y0 + dy
    plot(ptt, ptpup, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[1],
         lwd=lwd[1], ...)
    fillplot(ptt, ptplow, ptpup, col=cols[1],lwd=lwd[1])
    text(xlim[2]-eps, y, legend[ord[1]], adj=1, cex=cex[1])
    for (s in 2:S) {
      ptplow <- ptpup
      ptpup <- ptpup + ptp[1:nt,ord[s]]
      fillplot(ptt, ptplow, ptpup, col=cols[s],lwd=lwd[s])
      y0 <- y1
      dy <- ptp[nt, ord[s]]
      y <- y0 + 0.5*dy
      y1 <- y0 + dy
      text(xlim[2]-eps, y, legend[ord[s]], adj=1, cex=cex[s])
    }
    box()
  }
  else if (type=="separate") {
    if (missing(xlim)) xlim <- range(ptt)
    for (s in 1:S) {
      if (missing(ylim)) # each on its own y-scale
        plot(ptt, ptp[,s], type="s", xlim=xlim,  xlab=xlab, ylab=ylab, col=cols[s], lwd=lwd[s])
      else
        plot(ptt, ptp[,s], type="s", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=cols[s],
             lwd=lwd[s], ...)
      title(main=legend[s])
    }
  }
  return(invisible())
}
