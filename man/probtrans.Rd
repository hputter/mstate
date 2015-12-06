\name{probtrans}
\alias{probtrans}
\title{Compute subject-specific or overall transition probabilities with standard errors}
\description{
This function computes subject-specific or overall transition probabilities
in multi-state models. If requested, also standard errors are calculated.}
\usage{
probtrans(object, predt, direction=c("forward","fixedhorizon"),
          method=c("aalen","greenwood"), variance=TRUE, covariance=FALSE)
}
\arguments{
  \item{object}{\link{msfit} object containing estimated cumulative hazards
  for each of the transitions in the multi-state model and, if standard
  errors are requested, (co)variances of these cumulative hazards for each
  pair of transitions}
  \item{predt}{A positive number indicating the prediction time. This is
  either the time at which the prediction is made (if \code{direction}=
  \code{"forward"}) or the time for which the prediction is to be made (if
  \code{direction}=\code{"fixedhorizon"})}
  \item{direction}{One of \code{"forward"} (default) or \code{"fixedhorizon"},
  indicating whether prediction is forward or for a fixed horizon}
  \item{method}{A character string specifying the type of variances
  to be computed (so only needed if either \code{variance} or \code{covariance}
  is \code{TRUE}). Possible values are \code{"aalen"} or \code{"greenwood"}}
  \item{variance}{Logical value indicating whether standard errors are
  to be calculated (default is \code{TRUE})}
  \item{covariance}{Logical value indicating whether covariances of
  transition probabilities for different states are to be calculated
  (default is \code{FALSE})}
}
\details{
For details refer to de Wreede, Fiocco & Putter (2010).
}
\value{An object of class \code{"probtrans"}, which is a list of which item
  [[s]] contains a data frame with the estimated transition probabilities (and
  standard errors if \code{variance}=\code{TRUE}) from state s.
  If \code{covariance}=\code{TRUE}, item \code{varMatrix} contains an array
  of dimension K^2 x K^2 x (nt+1) (with K the number of states and nt the
  distinct transition time points); the time points correspond to those in
  the data frames with the estimated transition probabilities.
  Finally, there are items \code{trans}, \code{method}, \code{predt},
  \code{direction}, recording the transition matrix, and the method, predt
  and direction arguments used in the call to probtrans.
  Plot and summary methods have been defined for \code{"probtrans"} objects.
}
\references{
Andersen PK, Borgan O, Gill RD, Keiding N (1993). \emph{Statistical Models
Based on Counting Processes}. Springer, New York.

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
\author{Liesbeth de Wreede and Hein Putter \email{H.Putter@lumc.nl}}
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
HvH <- msfit(cx,newdata,trans=tmat)
# probtrans
pt <- probtrans(HvH,predt=0)
# predictions from state 1
pt[[1]]
}
\keyword{survival}
