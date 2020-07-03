\name{etm2msdata}
\alias{etm2msdata}
\alias{trans2tra}
\alias{msdata2etm}
\alias{tra2trans}
\title{Converts between etm and msdata format}
\description{
Converts multi-state data back and forth between etm and msdata formats.
Covariates have to be dealt with separately.
}
\usage{
msdata2etm(msdata, id)
trans2tra(trans)
tra2trans(tra)
etm2msdata(etmdata, id, tra)
}
\arguments{
  \item{msdata}{Multi-state data in \code{msdata} format, as used in \code{mstate}}
  \item{id}{Column name identifying the subject id}
  \item{trans}{Transition matrix in \code{mstate} format}
  \item{etmdata}{Multi-state data in \code{etm} format}
  \item{tra}{Transition matrix in \code{etm} format}
}
\details{
\code{msdata2etm} will convert from \code{msdata} format to \code{etm} format;
\code{etm2msdata} will convert from \code{etm} format to \code{msdata} format
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
# Transition matrix for illness-death model
tmat <- trans.illdeath()
# Data in wide format, for transition 1 this is dataset E1 of
# Therneau & Grambsch (T&G)
tg <- data.frame(id=1:6,illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
        dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
        x1=c(1,1,1,0,0,0),x2=c(6:1))
# Data in long format using msprep
tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
        data=tg,keep=c("x1","x2"),trans=tmat, id="id")
# Same thing in etm format
tra <- trans2tra(tmat)
tgetm <- msdata2etm(tglong, id="id")
# And back
etm2msdata(tgetm, id="id", tra=tra)
}
\keyword{datagen}
