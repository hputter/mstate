\name{paths}
\alias{paths}
\title{Find all possible trajectories through a given multi-state model}
\description{
For a given multi-state model, specified through a transition matrix,
\code{paths} recursively finds all the possible trajectories or paths
through that multi-state starting from a specified state. DO NOT USE
for reversible or cyclic multi-state models.}
\usage{paths(trans, start=1)}
\arguments{
  \item{trans}{The transition matrix describing the multi-state
  model, see \code{\link{msprep}}}
  \item{start}{The starting state for the trajectories}
}
\details{
The function is recursive. It starts in \code{start}, looks at what
states can be visited from \code{start}, and appends the results of
the next call to the current value (matrix). If the transition matrix
contains loops, the function will find infinitely many paths, so do not
use \code{paths} for reversible or cyclic multi-state models. A
warning is not yet incorporated!
}
\value{
A matrix, each row of which specifies a possible path through the
multi-state model.
}
\author{Hein Putter <H.Putter@lumc.nl>}
\examples{
tmat <- matrix(NA,5,5)
tmat[1,2:3] <- 1:2
tmat[1,5] <- 3
tmat[2,4:5] <- 4:5
tmat[3,4:5] <- 6:7
tmat[4,5] <- 8
paths(tmat)
paths(tmat, start=3)
}
\keyword{array}
