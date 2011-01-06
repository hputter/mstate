\name{Cuminc}
\alias{Cuminc}
\title{Calculate nonparametric cumulative incidence functions and
associated standard errors}
\description{
This function computes nonparametric cumulative incidence functions and
associated standard errors for each value of a group variable.}
\usage{
Cuminc(time, status, data, group, failcodes, na.status=c("remove","extra"),
        variance=TRUE)
}
\arguments{
  \item{time}{Either 1) a numeric vector containing the failure times
  or 2) a string containing the column name indicating these failure times}
  \item{status}{Either 1) a numeric, factor or character vector containing
  the failure codes or 2) a string containing the column name indicating
  these failure codes}
  \item{data}{When appropriate, a data frame containing \code{time},
  \code{status} and/or \code{group} variables}
  \item{group}{Optionally, name of column in data indicating a
  grouping variable; cumulative incidence functions are calculated
  for each value or level of \code{group}. If missing no groups are
  considered}
  \item{failcodes}{A vector indicating which values of \code{status} are
  considered as different causes of failure; other values of \code{status}
  are considered as censorings. If missing and \code{status} is numeric,
  it is assumed that 0 is censoring and all other values indicate failcodes;
  if missing and \code{status} is character or factor, then it is assumed
  that each of the levels/values of \code{status} is a cause of failure}
  \item{na.status}{One of \code{"remove"} (default) or \code{"extra"},
  indicating whether subjects with missing cause of failure should be
  removed or whether missing cause of failure should be treated as a
  separate cause of failure}
  \item{variance}{Logical value, indicating whether the standard errors
  of the cumulative incidences should be output (\code{TRUE}, the default)
  or not}
}
\details{
The estimated cumulative incidences are as described in Putter,
Fiocco & Geskus (2007); the standard errors are the square roots of
the Greenwood variance estimators, see eg. Andersen, Borgan, Gill & Keiding
(1993), de Wreede, Fiocco & Putter (2009), and they correspond to the
variances in eg. Marubini & Valsecchi (1997). In case of no
censoring, the estimated cumulative incidences and variances reduce
to simple binomial frequencies and their variances.
}
\value{A data frame containing the estimated failure-free
probabilities and cumulative incidences and their standard errors.
The names of the dataframe are \code{time}, \code{Surv},
\code{seSurv}, and \code{cuminc} and \code{secuminc} followed by
the values or levels of the \code{failcodes}. If \code{group} was
specified, a \code{group} variable is included with the same name and
values/levels as the original grouping variable, and with estimated
cumulative incidences (SE) for each value/level of \code{group}.
}
\references{
Andersen PK, Borgan O, Gill RD, Keiding N (1993). \emph{Statistical Models
Based on Counting Processes}. Springer, New York.

Marubini E, Valsecchi MG (1995). \emph{Analysing Survival Data from
Clinical Trials and Observational Studies}. Wiley, New York.

Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
2389--2430.

de Wreede L, Fiocco M, Putter H (2009). The mstate package for
estimation and prediction in non- and semi-parametric multi-state
models. Submitted. \url{www.msbi.nl/multistate}.
}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
### These data were used in Putter, Fiocco & Geskus (2007)
data(aidssi)
ci <- Cuminc(time=aidssi$time, status=aidssi$status)
head(ci); tail(ci)
ci <- Cuminc(time="time", status="status", data=aidssi, group="ccr5")
head(ci); tail(ci)

### Some fake data
fake <- data.frame(surv=c(seq(2,10,by=2),seq(1,13,by=3),seq(1,9,by=2),seq(1,13,by=3)),
                    stat=rep(0:3,5),Tstage=c(1:4,rep(1:4,rep(4,4))))
fake$stat[fake$stat==0 & fake$Tstage==2] <- 3
fake$stat[fake$stat==3 & fake$Tstage==1] <- 2
fake
Cuminc(time="surv", status="stat", data=fake)
# If we remove all entries with status=0,
# we should get binomial sample probabilities and corresponding SEs
fake0 <- fake[fake$stat!=0,]
Cuminc(time="surv", status="stat", data=fake0)
# Use failcodes
Cuminc(time="surv", status="stat", data=fake, failcodes=c(1,3))
# Make grouping variable and status variable a factor
fake$Tstage <- factor(fake$Tstage,labels=c("T1","T2","T3","T4"))
fake$stat <- factor(fake$stat,levels=0:3,
            labels=c("eventfree","event1","event2","event3"))
Cuminc(time="surv", status="stat", data=fake, group="Tstage")
# (The warnings are a result of the fact that some failure causes
# do not occur for some values/levels of the grouping variable)
# The previous command didn't do what we wanted because we didn't
# tell Cuminc that eventfree is not an event
Cuminc(time="surv", status="stat", data=fake, group="Tstage",
        failcodes=c("event1","event2","event3"))
}
\keyword{survival}
