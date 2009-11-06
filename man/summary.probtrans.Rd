\name{summary.probtrans}
\alias{summary.probtrans}
\title{Summary method for a probtrans object}
\description{
Summary method for an object of class 'probtrans'. It prints a selection
of the estimated transition probabilities, and, if requested, also of the
variances.}
\usage{
\S3method{summary}{probtrans}(object,from,complete=FALSE,variance=TRUE,...)
}
\arguments{
  \item{object}{Object of class 'probtrans', containing estimated
  transition probabilities from and to all states in a multi-state model}
  \item{from}{Specifies from which state the transition probabilities
  are to be printed. Should be subset of 1:S, with S the number of states
  in the multi-state model. If missing, transition probabilities are
  printed from all starting states}
  \item{complete}{Whether or not the complete estimated transition
  probabilities should be printed (\code{TRUE}) or not (\code{FALSE});
  default is \code{FALSE}, in which case the estimated transition
  probilities will be printed for the first and last 6 time points of
  each starting state (or all when there are at most 12 of these time
  points}
  \item{variance}{Whether or not the standard errors of the estimated
  transition probabilities should be printed; default is \code{TRUE}}
  \item{\dots}{Further arguments to summary}
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
# default summary
summary(pt)
# summary without standard errors
summary(pt,variance=FALSE)
}
\keyword{print}
