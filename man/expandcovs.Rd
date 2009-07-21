\name{expand.covs}
\alias{expand.covs}
\title{Expand covariates in dataset in long format}
\description{
Given a multi-state dataset in long format, and one or more
covariates, this function adds transition-specific covariates,
expanding the original covariate(s), to the dataset. The original
dataset with the transition-specific covariates appended is returned.
}
\usage{expand.covs(data,trans,covs,append=TRUE,longnames=TRUE)}
\arguments{
  \item{data}{Data frame in long format, such as output by
  \code{\link{msprep}}}
  \item{trans}{Transition matrix describing the states and
  transitions in the multi-state model. See \code{trans} in
  \code{\link{msprep}} for more detailed information}
  \item{covs}{A character vector containing the names of
  the covariates in \code{data} to be expanded}
  \item{append}{Logical value indicating whether or not the design
  matrix for the expanded covariates should be appended to the data
  (default=\code{TRUE})}
  \item{longnames}{Logical value indicating whether or not the labels
  are to be used for the names of the expanded covariates that are
  categorical (default=\code{TRUE}); in case of \code{FALSE} numbers
  from 1 up to the number of contrasts are used}
}
\details{
For a given basic covariate \code{Z}, the transition-specific
covariate for transition \code{s} is called \code{Z.s}. The concept of
transition-specific covariates in the context of multi-state models
was introduced by Andersen, Hansen & Keiding (1991), see also Putter, 
Fiocco & Geskus (2007). It is
only unambiguously defined for numeric covariates or for explicit
codings. Then it will take the value 0 for all rows in the long
format dataframe for which \code{trans} does not equal \code{s}.
For the rows for which \code{trans} equals \code{s}, the original
value of \code{Z} is copied. In \code{expand.covs}, when a given
covariate is a factor, it will be expanded on the design matrix
given by \code{\link[stats:model.matrix]{model.matrix}}. Missing values
in the basic covariates are allowed and result in missing values
in the expanded covariates.
}
\value{
A data frame containing the design matrix for the transition-
specific covariates, either on its own (\code{append}=\code{FALSE})
or appended to the data (\code{append}=\code{TRUE}).
}
\references{
Andersen PK, Hansen LS, Keiding N (1991). Non- and semi-parametric
estimation of transition probabilities from censored observation of a
non-homogeneous Markov process. \emph{Scandinavian Journal of Statistics}
\bold{18}, 153--167.

Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
2389--2430.
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
# transition matrix for illness-death model
tmat <- trans.illdeath()
# small data set in wide format
tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
        dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
        x1=c(1,1,1,2,2,2),x2=c(6:1))
tg$x1 <- factor(tg$x1,labels=c("male","female"))
# data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),
        status=c(NA,"ills","ds"),data=tg,
        keep=c("x1","x2"),trans=tmat)
# expanded covariates
expand.covs(tglong,tmat,c("x1","x2"),append=FALSE)
expand.covs(tglong,tmat,"x1")
expand.covs(tglong,tmat,"x1",longnames=FALSE)
}
\keyword{datagen}
