\name{transMat}
\alias{transMat}
\alias{trans.illdeath}
\alias{trans.comprisk}
\title{Define transition matrix for multi-state model}
\description{
Define transition matrices for multi-state model. Specific functions for
defining such transition matrices are pre-defined for common multi-state models
like the competing risks model and the illness-death model.}
\usage{
transMat(x, names)
trans.comprisk(K, names)
trans.illdeath(names)
}
\arguments{
  \item{x}{List of possible transitions; x[[i]] consists of a vector of
  state numbers reachable from state i}
  \item{K}{The number of competing risks}
  \item{names}{A character vector containing the names of either the
  competing risks or the states in the multi-state model specified by
  the competing risks or illness-death model. \code{names} should have
  the same length as the list \code{x} (for \code{transMat}), or either
  \code{K} or \code{K}+1 (for \code{trans.comprisk}), or 3 (for
  \code{trans.illdeath})}
}
\details{
If \code{names} is missing, the names \code{"eventfree"}, \code{"cause1"},
etcetera are assigned in \code{trans.comprisk}, or \code{"healthy"},
\code{"illness"}, \code{"death"} in \code{trans.illdeath}.}
\value{
A transition matrix describing the states and transitions in the
multi-state model.
}
\author{Steven McKinney <smckinney@bccrc.ca>; Hein Putter <H.Putter@lumc.nl>}
\examples{
transMat(list(c(2, 3), c(), c(1, 2)),
	names = c("Disease-free", "Death", "Relapsed"))
tmat <- transMat(x = list( c(2, 3), c(1, 3), c() ),
                 names = c("Normal", "Low", "Death"))
tmat
transListn <- list("Normal" = c(2, 3), "Low" = c(1, 3), "Death" = c())
transMat(transListn)
trans.comprisk(3)
trans.comprisk(3,c("c1","c2","c3"))
trans.comprisk(3,c("nothing","c1","c2","c3"))
trans.illdeath()
trans.illdeath(c("nothing","ill","death"))
}
\keyword{array}
