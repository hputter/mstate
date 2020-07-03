\name{transhelp}
\alias{to.trans2}
\alias{trans2Q}
\alias{absorbing}
\alias{is.circular}
\title{Help functions for transition matrix}
\description{
Help functions to get insight into the structure of a transition matrix.}
\usage{
to.trans2(trans)
trans2Q(trans)
absorbing(trans)
is.circular(trans)
}
\arguments{
  \item{trans}{Transition matrix, for instance produced by \code{transMat}),
  \code{trans.comprisk}, or \code{trans.illdeath}}
}
\details{
Function \code{to.trans2} simply lists the transitions in \code{trans} in a data frame; function
\code{trans2Q} converts \code{trans} to a \code{Q} matrix, the (j,k)th element of which contains
the (shortest) number of transitions needed to travel from the jth to the kth state; function
\code{absorbing} returns a vector (named if \code{trans} contains row or columnc names) with the
state numbers that are absorbing; function \code{is.circular} returns (a Boolean) whether the
transition matrix specified in \code{trans} is circular or not.}
\value{
See details.
}
\author{Hein Putter <H.Putter@lumc.nl>}
\examples{
# Irreversible illness-death model
tmat <- trans.illdeath(c("Healthy", "Illness", "Death"))
tmat
to.trans2(tmat)
trans2Q(tmat)
absorbing(tmat)
is.circular(tmat)
# Reversible illness-death model
tmat <- transMat(x = list( c(2, 3), c(1, 3), c() ),
                 names = c("Healthy", "Illness", "Death"))
tmat
to.trans2(tmat)
trans2Q(tmat)
absorbing(tmat)
is.circular(tmat)
}
\keyword{univar}
