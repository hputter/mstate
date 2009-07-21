\name{trans.illdeath}
\alias{trans.illdeath}
\alias{trans.comprisk}
\title{Define transition matrix for competing risks or illness-death model}
\description{
Define transition matrices for common multi-state models like the
competing risks model and the illness-death model.}
\usage{
trans.comprisk(K, names)
trans.illdeath(names)
}
\arguments{
  \item{K}{The number of competing risks}
  \item{names}{A character vector containing the names of either the
  competing risks or the states in the multi-state model specified by
  the competing risks or illness-death model. The length of \code{names}
  should be either \code{K} or \code{K}+1 (for \code{trans.comprisk}),
  or 3 (for \code{trans.illdeath})}
}
\details{
If \code{names} is missing, the names \code{"eventfree"}, \code{"cause1"},
etcetera are assigned in \code{trans.comprisk}, or \code{"healthy"},
\code{"illness"}, \code{"death"} in \code{trans.illdeath}.}
\value{
A transition matrix describing the states and transitions in the
multi-state model.
}
\author{Hein Putter <H.Putter@lumc.nl>}
\examples{
trans.comprisk(3)
trans.comprisk(3,c("c1","c2","c3"))
trans.comprisk(3,c("nothing","c1","c2","c3"))
trans.illdeath()
trans.illdeath(c("nothing","ill","death"))
}
\keyword{array}
