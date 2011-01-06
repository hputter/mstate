\name{msprep}
\alias{msprep}
\title{Function to prepare dataset for multi-state modeling in long
format from dataset in wide format}
\description{
This function converts a dataset which is in wide format (one subject
per line, multiple columns indicating time and status for different
states) into a dataset in long format (one line for each transition
for which a subject is at risk). Selected covariates are replicated
per subjects.}
\usage{
msprep(time, status, data, trans, start, id, keep)
}
\arguments{
  \item{time}{Either 1) a matrix or data frame of dimension n x S
  (n being the number of individuals and S the number of states in
  the multi-state model), containing the times
  at which the states are visited or last follow-up time, or
  2) a character vector of length S containing the column
  names indicating these times. In the latter cases, some elements
  of \code{time} may be NA, see Details}
  \item{status}{Either 1) a matrix or data frame  of dimension n x S,
  containing, for each
  of the states, event indicators taking the value 1 if the state is
  visited or 0 if it is  not (censored), or
  2) a character vector of length S containing the column
  names indicating these status variables. In the latter cases,
  some elements of \code{status} may be NA, see Details}
  \item{data}{Data frame in wide format in which to interpret
  \code{time}, \code{status}, \code{id} or \code{keep}, if appropriate}
  \item{trans}{Transition matrix describing the states and transitions
  in the multi-state model. If S is the number of states in the
  multi-state model, \code{trans} should be an S x S matrix,
  with (i,j)-element a positive integer if a transition from i to j
  is possible in the multi-state model, \code{NA} otherwise. In particular,
  all diagonal elements should be \code{NA}. The
  integers indicating the possible transitions in the multi-state
  model should be sequentially numbered, 1,...,K, with K the number
  of transitions}
  \item{start}{List with elements \code{state} and \code{time},
  containing starting states and times of the subjects in the data.
  Default is \code{NULL}, in which case all subjects start in state 1
  at time 0. If a single state and time are given this state and time
  is used for all subjects, otherwise the length of \code{state} and
  \code{time} should equal the number of subjects in \code{data}}
  \item{id}{Either 1) a vector of length n containing the subject identifications,
  or 2) a character string indicating the column name containing
  these subject ids. If not provided, \code{"id"} will be assigned
  with values 1,...,n}
  \item{keep}{Either 1) a data frame or matrix with n rows or a numeric
  or factor vector of length n
  containing covariate(s) that need to be retained in the output
  dataset, or 2) a character vector containing the column names of these
  covariates in \code{data}}
}
\details{
For \code{msprep}, the transition matrix should correspond to an
irreversible acyclic Markov chain. In particular, on the diagonals
only \code{NA}s are allowed.

The transition matrix, if irreversible and acyclic, will have starting
states, i.e. states into which no transitions are possible. For these
starting states \code{NA}s
are allowed in the \code{time} and \code{status} arguments, either as
columns, when specified as matrix or data frame, or as elements of the
character vector when specified as character vector.

The function \code{msprep} uses a recursive algorithm through calls
to the recursive function \code{msprepEngine}. First, with the
current transition matrix, all starting
states are detected (defined as states into which there are no
transitions). For each of these starting states, all subjects starting from
that state are selected and for each subject the next visited state
is detected by looking at all transitions from that starting state
and determining the smallest transition time with \code{status}=1.
The recursive \code{msprepEngine} is called again with the starting
states deleted from the transition matrix and with subjects deleted
that either reached an absorbing state or that were censored. For
the remaining subjects the starting states and times are updated in
the next call. Datasets returned from the \code{msprepEngine} calls
are appended to the current dataset in long format and finally sorted.

A warning is issued for a subject, if multiple transitions exist with
the same smallest transition time (and \code{status}=0). In such
cases the next transition cannot be determined unambiguously, and the
state with the smallest number is chosen. In our experience,
occasionally the shortest transition time has \code{status}=0, while
a higher transition time has \code{status}=1. Then this larger
transition time and the corresponding transition is selected. No
warning is issued for these data inconsistencies. 
}
\value{An object of class \code{"msdata"}, which is a data frame in long
(counting process) format containing the subject id, the covariates
(replicated per subject), and
  \item{from}{the starting state}
  \item{to}{the receiving state}
  \item{trans}{the transition number}
  \item{Tstart}{the starting time of the transition}
  \item{Tstop}{the stopping time of the transition}
  \item{status}{status variable, with 1 indicating an
  event (transition), 0 a censoring}
  The \code{"msdata"} object has the transition matrix as \code{"trans"}
  attribute.
  }
\references{
Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine}
\bold{26}, 2389--2430.
}
\author{Hein Putter \email{H.Putter@lumc.nl} and Marta Fiocco}
\examples{
# transition matrix for illness-death model
tmat <- trans.illdeath()
# some data in wide format
tg <- data.frame(stt=rep(0,6),sts=rep(0,6),
        illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
        dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
        x1=c(1,1,1,2,2,2),x2=c(6:1))
tg$x1 <- factor(tg$x1,labels=c("male","female"))
tg$patid <- factor(2:7,levels=1:8,labels=as.character(1:8))
# define time, status and covariates also as matrices
tt <- matrix(c(rep(NA,6),tg$illt,tg$dt),6,3)
st <- matrix(c(rep(NA,6),tg$ills,tg$ds),6,3)
keepmat <- data.frame(gender=tg$x1,age=tg$x2)
# data in long format using msprep
msprep(time=tt,status=st,trans=tmat,keep=as.matrix(keepmat))
msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),data=tg,
		id="patid",keep=c("x1","x2"),trans=tmat)
# Patient no 5, 6 now start in state 2 at time t=4 and t=10
msprep(time=tt,status=st,trans=tmat,keep=keepmat,
        start=list(state=c(1,1,1,1,2,2),time=c(0,0,0,0,4,10)))
}
\keyword{datagen}
