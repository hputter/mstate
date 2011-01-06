\name{summary.msfit}
\alias{summary.msfit}
\title{Summary method for an msfit object}
\description{
Summary method for an object of class 'msfit'. It prints a selection
of the estimated cumulative transition intensities, and, if requested,
also of the (co)variances.}
\usage{
\S3method{summary}{msfit}(object,complete=FALSE,variance=FALSE,...)
}
\arguments{
  \item{object}{Object of class 'msfit', containing estimated cumulative
  transition intensities for all transitions in a multi-state model}
  \item{complete}{Whether or not the complete estimated cumulative
  transition intensities should be printed (\code{TRUE}) or not
  (\code{FALSE}); default is \code{FALSE}, in which case the estimated
  cumulative transition hazards will be printed for the first and last
  6 time points of each transition (or all when there are at most 12
  of these time points}
  \item{variance}{Whether or not the (co)variances of the estimated
  cumulative transition intensities should be printed; default is
  \code{FALSE}}
  \item{\dots}{Further arguments to summary}
}
\value{No return value}
\seealso{\code{\link{msfit}}}
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
print(msf)
# standard summary
summary(msf)
# including variances and covariances
summary(msf,variance=TRUE)
}
\keyword{print}
