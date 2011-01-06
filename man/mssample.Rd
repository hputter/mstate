\name{mssample}
\alias{mssample}
\title{Sample paths through a multi-state model}
\description{
Given cumulative transition hazards sample paths through the multi-state model.
}
\usage{
mssample(Haz, trans, history=list(state=1,time=0,tstate=NULL),
    beta.state=NULL, clock=c("forward","reset"),
    output=c("state","path","data"),
    tvec, cens=NULL, M=10, do.trace=NULL)
}
\arguments{
  \item{Haz}{Cumulative hazards to be sampled from. These should be
  given as a data frame with columns \code{time}, \code{Haz}, \code{trans},
  for instance as the \code{Haz} list element given by \code{\link{msfit}}.}
  \item{trans}{Transition matrix describing the multi-state model.  See
  \code{trans} in \code{\link{msprep}} for more detailed information}
  \item{history}{A list with elements \code{state}, specifying the starting
  state(s), \code{time}, the starting time(s), and \code{tstate}, a numeric
  vector of length the number of states, specifying at what times states
  have been visited, if appropriate. The default of \code{tstate} is
  \code{NULL}; more information can be found under Details.
  
  The elements \code{state} and \code{time} may either be scalars or vectors,
  in which case different sampled paths may start from different states or
  at different times. By default, all sampled paths start from state 1 at
  time 0.}
  \item{beta.state}{A matrix of dimension (no states) x (no transitions)
  specifying estimated effects of times at which earlier states were reached
  on subsequent transitions. If these are not in the model, the value
  \code{NULL} (default) suffices; more information can be found under Details}
  \item{clock}{Character argument, either \code{"forward"} (default) or
  \code{"reset"}, specifying whether the time-scale of the cumulative hazards
  is in forward time (\code{"forward"}) or duration in the present state
  (\code{"reset"})}
  \item{output}{One of \code{"state"}, \code{"path"}, or \code{"data"},
  specifying whether states, paths, or data should be output.}
  \item{tvec}{A numeric vector of time points at which the states or paths
  should be evaluated. Ignored if \code{output}=\code{"data"}}
  \item{cens}{An independent censoring distribution, given as a data frame
  with time and Haz}
  \item{M}{The number of sampled trajectories through the multi-state model.
  The default is 10, since the procedure can become quite time-consuming}
  \item{do.trace}{An integer, specifying that the replication number should be
  written to the console every \code{do.trace} replications. Default is
  \code{NULL} in which case no output is written to the console during the
  simulation}
}
\details{
The procedure is described in detail in Fiocco, Putter & van Houwelingen (2008).
The argument \code{beta.state} and the element \code{tstate} from the
argument \code{history} are meant to incorporate situations where the time at
which some previous states were visited may affect future transition rates.
The relation between time of visit of state \code{s} and transition \code{k}
is assumed to be linear on the log-hazards; the corresponding regression
coefficient is to be supplied as the (s,k)-element of \code{beta.state},
which is 0 if no such effect has been included in the model. If no such effects
are present, then \code{beta.state}=\code{NULL} (default) suffices. In the
\code{tstate} element of \code{history}, the \code{s}-th element is the time
at which state \code{s} was visited. This is only relevant for states which
have been visited prior to the beginning of sampling, i.e. before the
\code{time} element of \code{history}; the elements of \code{tstate} are
internally updated when in the sampling process new states are visited
(only if \code{beta.state} is not \code{NULL} to avoid unnecessary
computations).
}
\value{
M simulated paths through the multi-state model given by \code{trans} and
\code{Haz}. It is either a data frame with columns \code{time}, \code{pstate1},
..., \code{pstateS} for S states when \code{output="state"}, or with columns
\code{time}, \code{ppath1},..., \code{ppathP} for the P paths specified in
\code{\link{paths}}(trans) when \code{output="path"}. When \code{output="data"},
the sampled paths are stored in an \code{"msdata"} object, a data frame in long
format such as that obtained by \code{\link{msprep}}. This may be useful for
(semi-)parametric bootstrap procedures, in which case \code{cens} may be used
as censoring distribution (assumed to be independent of all transition times
and independent of any covariates).
}
\references{
Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
proportional hazards regression and simulation-based prediction for
multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.
}
\author{Marta Fiocco, Hein Putter \email{H.Putter@lumc.nl}}
\examples{
# transition matrix for illness-death model
tmat <- trans.illdeath()
# data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (T&G)
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
        dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
        x1=c(1,1,1,0,0,0),x2=c(6:1))
# data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
        data=tg,keep=c("x1","x2"),trans=tmat)
# expanded covariates
tglong <- expand.covs(tglong,c("x1","x2"))
# Cox model with different covariate
cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
        data=tglong,method="breslow")
# new data, to check whether results are the same for transition 1 as T&G
newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
fit <- msfit(cx,newdata,trans=tmat)
tv <- unique(fit$Haz$time)
# mssample
set.seed(1234)
mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100)
set.seed(1234)
paths(tmat)
mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100,output="path")
set.seed(1234)
mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100,output="data",do.trace=25)
}
\keyword{datagen}
