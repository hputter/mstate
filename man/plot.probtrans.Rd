\name{plot.probtrans}
\alias{plot.probtrans}
\title{Plot method for a probtrans object}
\description{
Plot method for an object of class 'probtrans'. It plots the transition
probabilities as estimated by \code{\link{probtrans}}.}
\usage{
\S3method{plot}{probtrans}(x,from=1,
    type=c("stacked","filled","single","separate"),ord,
    cols,xlab="Time",ylab="Probability",ylim,lwd,lty,cex,
    legend,legend.pos,bty="o",...)
}
\arguments{
  \item{x}{Object of class 'probtrans', containing estimated transition
  probabilities}
  \item{from}{The starting state from which the probabilities are used to
  plot}
  \item{type}{One of \code{"stacked"} (default), \code{"filled"},
  \code{"single"} or \code{"separate"}; in case of \code{"stacked"}, the
  transition probabilities are stacked and the distance between two
  adjacent curves indicates the probability, this is also true for
  \code{"filled"}, but the space between adjacent curves are filled,
  in case of \code{"single"}, the probabilities are shown as different
  curves in a single plot, in case of \code{"separate"}, separate plots
  are shown for the estimated transition probabilities}
  \item{ord}{A vector of length equal to the number of states, specifying
  the order of plotting in case type=\code{"stacked"} or \code{"filled"}}
  \item{cols}{A vector specifying colors for the different transitions;
  default is 1:K (K no of transitions), when type=\code{"single"}, and
  1 (black), when type=\code{"separate"}}
  \item{xlab}{A title for the x-axis; default is \code{"Time"}}
  \item{ylab}{A title for the y-axis; default is
  \code{"Probability"}}
  \item{ylim}{The y limits of the plot(s); if ylim is specified for
  type="separate", then all plots use the same ylim for y limits}
  \item{lwd}{The line width, see \code{\link{par}}; default is 1}
  \item{lty}{The line type, see \code{\link{par}}; default is 1}
  \item{cex}{}
  \item{legend}{Character vector of length equal to the number of
  transitions, to be used in a legend; if missing, numbers will be used;
  this and the legend arguments following are ignored when type="separate"}
  \item{legend.pos}{The position of the legend, see \code{\link{legend}};
  default is \code{"topleft"}}
  \item{bty}{The box type of the legend, see \code{\link{legend}}}
  \item{\dots}{Further arguments to plot}
}
\value{No return value}
\seealso{\code{\link{probtrans}}}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
# transition matrix for illness-death model
tmat <- trans.illdeath()
# data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (2000)
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
        dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
        x1=c(1,1,1,0,0,0),x2=c(6:1))
# data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
		data=tg,keep=c("x1","x2"),trans=tmat)
# events
events(tglong)
table(tglong$status,tglong$to,tglong$from)
# expanded covariates
tglong <- expand.covs(tglong,c("x1","x2"))
# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
	data=tglong,method="breslow")
summary(cx)
# new data, to check whether results are the same for transition 1 as
# those in appendix E.1 of Therneau & Grambsch (2000)
newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
msf <- msfit(cx,newdata,trans=tmat)
# probtrans
pt <- probtrans(msf,predt=0)
# default plot
plot(pt,ord=c(2,3,1),lwd=2,cex=0.75)
# filled plot
plot(pt,type="filled",ord=c(2,3,1),lwd=2,cex=0.75)
# single plot
plot(pt,type="single",lwd=2,col=rep(1,3),lty=1:3,legend.pos=c(8,1))
# separate plots
par(mfrow=c(2,2))
plot(pt,type="sep",lwd=2)
par(mfrow=c(1,1))
}
\keyword{hplot}
