\name{msboot}
\alias{msboot}
\title{Bootstrap function in multi-state models}
\description{
A generic nonparametric bootstrapping function for multi-state models.
}
\usage{
msboot(theta, data, trans, B=5, id="id", verbose=0, ...)
}
\arguments{
  \item{theta}{A function of \code{data}, \code{trans}, and perhaps other
  arguments, returning the value of the statistic to be bootstrapped; the output
  of theta should be a scalar or numeric vector}
  \item{data}{The original data in long format, such as output from
  \code{\link{msprep}}}
  \item{trans}{Transition matrix describing the states and transitions
  in the multi-state model. See \code{trans} in \code{\link{msprep}} for
  more detailed information}
  \item{B}{The number of bootstrap replications; the default is taken
  to be quite small (5) since bootstrapping can be time-consuming}
  \item{id}{Character string indicating which column identifies the
  subjects to be resampled}
  \item{verbose}{The level of output; default 0 = no output,
  1 = print the replication}
  \item{...}{Any further arguments to the function \code{theta}}
}
\details{
The function \code{msboot} samples randomly with replacement subjects
from the original dataset \code{data}. The individuals are identified
with \code{id}, and bootstrap datasets are produced by concatenating
all selected rows. 
}
\value{
Matrix of dimension (length of output of theta) x B, with b'th
column being the value of theta for the b'th bootstrap dataset
}
\references{
Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
proportional hazards regression and simulation-based prediction for
multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.
}
\author{Marta Fiocco, Hein Putter <H.Putter@lumc.nl>}
\examples{
tmat <- trans.illdeath()
data(ebmt1)
covs <- c("score","yrel")
msebmt <- msprep(time=c(NA,"rel","srv"),status=c(NA,"relstat","srvstat"),
		data=ebmt1,id="patid",keep=covs,trans=tmat)
# define a function (this one returns vector of regression coef's)
regcoefvec <- function(data,trans) {
  cx <- coxph(Surv(Tstart,Tstop,status)~score+strata(trans),
          data=data,method="breslow")
  return(coef(cx))
}
regcoefvec(msebmt,tmat)
set.seed(1234)
msboot(theta=regcoefvec,data=msebmt,trans=tmat,id="patid")
}
\keyword{datagen}
