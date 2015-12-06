\name{msfit}
\alias{msfit}
\title{Compute subject-specific transition hazards with (co-)variances}
\description{
This function computes subject-specific or overall cumulative transition
hazards for each of the possible transitions in the multi-state model.
If requested, also the variances and covariances of the estimated cumulative
transition hazards are calculated.}
\usage{
msfit(object, newdata, variance=TRUE, vartype=c("aalen","greenwood"), trans)
}
\arguments{
  \item{object}{A \code{\link[survival:coxph]{coxph}} object
  describing the fit of the multi-state model}
  \item{newdata}{A data frame with the same variable names as those
  that appear in the \code{coxph} formula. Its use is somewhat
  different from \code{\link[survival:survfit]{survfit}}. See Details.
  The argument \code{newdata} may be omitted only if the right hand
  side of the formula in the \code{coxph} object is \code{~strata(trans)}}
  \item{variance}{A logical value indicating whether the
  (co-)variances of the subject-specific transition hazards should be
  computed. Default is \code{TRUE}}
  \item{vartype}{A character string specifying the type of variances
  to be computed (so only needed if \code{variance}=\code{TRUE}).
  Possible values are \code{"aalen"} or \code{"greenwood"}}
  \item{trans}{Transition matrix describing the states and transitions
  in the multi-state model. See \code{trans} in \code{\link{msprep}} for
  more detailed information}
}
\details{
The data frame needs to have one row for each transition in the
multi-state model. An additional column \code{strata} (numeric) is
needed to describe for each transition to which stratum it belongs.
The name has to be \code{strata}, even if in the original \code{coxph}
call another variable was used.
For details refer to de Wreede, Fiocco & Putter (2010). So far, the
results have been checked only for the \code{"breslow"} method of
dealing with ties in \code{\link[survival:coxph]{coxph}}, so this
is recommended.
}
\value{An object of class \code{"msfit"}, which is a list containing
  \item{Haz }{A data frame with \code{time}, \code{Haz},
  \code{trans}, containing the estimated subject-specific hazards for
  each of the transitions in the multi-state model}
  \item{varHaz }{A data frame
  with \code{time}, \code{Haz}, \code{trans1}, \code{trans2}
  containing the variances (\code{trans1}=\code{trans2}) and
  covariances (\code{trans1}<\code{trans2}) of the estimated
  hazards. This element is only returned when
  \code{variance}=\code{TRUE}}
  \item{trans}{The transition matrix used}
}
\references{
Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine}
\bold{26}, 2389--2430.

Therneau TM, Grambsch PM (2000). \emph{Modeling Survival Data:
Extending the Cox Model}. Springer, New York.

de Wreede LC, Fiocco M, and Putter H (2010). The mstate package
for estimation and prediction in non- and semi-parametric multi-state
and competing risks models. \emph{Computer Methods and Programs in
Biomedicine} \bold{99}, 261--274.

de Wreede LC, Fiocco M, and Putter H (2011). mstate: An R Package
for the Analysis of Competing Risks and Multi-State Models. \emph{Journal of
Statistical Software}, Volume 38, Issue 7.
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\seealso{\code{\link{plot.msfit}}}
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
msfit(cx,newdata,trans=tmat)
}
\keyword{survival}
