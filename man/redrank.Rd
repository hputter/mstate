\name{redrank}
\alias{redrank}
\title{Reduced rank proportional hazards model for competing risks and
multi-state models}
\description{
This function estimates regression coefficients in reduced rank
proportional hazards models for competing risks and multi-state models.}
\usage{
redrank(redrank, full = ~1, data, R, strata = NULL, Gamma.start,
            method = "breslow", eps = 1e-5, print.level = 1)
}
\arguments{
  \item{redrank}{Survival formula, starting with either Surv(time,status) ~
  or with Surv(Tstart,Tstop,status) ~, followed by a formula containing
  covariates for which a reduced rank estimate is to be found}
  \item{full}{Optional, formula specifying that part which needs to be
  retained in the model (so not subject to reduced rank)}
  \item{data}{Object of class 'msdata', as prepared for instance by
  \code{\link{msprep}}, in which to interpret the \code{redrank} and,
  optionally, the \code{full} formulas}
  \item{R}{Numeric, indicating the rank of the solution}
  \item{strata}{Name of covariate to be used inside the
  \code{\link[survival:strata]{strata}} part of
  \code{\link[survival:coxph]{coxph}}}
  \item{Gamma.start}{A matrix of dimension K x R, with K the number of
  transitions and R the rank, to be used as starting value}
  \item{method}{The method for handling ties in
  \code{\link[survival:coxph]{coxph}}}
  \item{eps}{Numeric value determining when the iterative algorithm
  is finished (when for two subsequent iterations the difference in
  log-likelihood is smaller than \code{eps})}
  \item{print.level}{Determines how much output is written to the
  screen; 0: no output, 1: iterations, for each iteration solutions of
  Alpha, Gamma, log-likelihood, 2: also the Cox regression results}
}
\details{
For details refer to Fiocco, Putter & van Houwelingen (2005, 2008).
}
\value{A list with elements
  \item{Alpha}{the Alpha matrix}
  \item{Gamma}{the Gamma matrix}
  \item{Beta}{the Beta matrix corresponding to \code{covariates}}
  \item{Beta2}{the Beta matrix corresponding to \code{fullcovs}}
  \item{cox.itr1}{the \code{\link[survival:coxph]{coxph}} object
  resulting from the last call giving \code{Alpha}}
  \item{alphaX}{the matrix of prognostic scores given by
  \code{Alpha}, n x R, with n number of subjects}
  \item{niter}{the number of iterations needed to reach convergence}
  \item{df}{the number of regression parameters estimated}
  \item{loglik}{the log-likelihood}
}
\references{
Fiocco M, Putter H, van Houwelingen JC (2005). Reduced rank
proportional hazards model for competing risks. \emph{Biostatistics}
\bold{6}, 465--478.

Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
proportional hazards regression and simulation-based prediction for
multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.

Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine}
\bold{26}, 2389--2430.
}
\author{Marta Fiocco and Hein Putter \email{H.Putter@lumc.nl}}
\examples{
# Based on Fiocco, Putter & van Houwelingen (2005)
data(ebmt2)
# transition matrix for competing risks
tmat <- trans.comprisk(6,names=c("Relapse","GvHD","Bacterial","Viral","Fungal","Other"))
# preparing long dataset
ebmt2$stat1 <- as.numeric(ebmt2$status==1)
ebmt2$stat2 <- as.numeric(ebmt2$status==2)
ebmt2$stat3 <- as.numeric(ebmt2$status==3)
ebmt2$stat4 <- as.numeric(ebmt2$status==4)
ebmt2$stat5 <- as.numeric(ebmt2$status==5)
ebmt2$stat6 <- as.numeric(ebmt2$status==6)
covs <- c("dissub","match","tcd","year","age")
ebmtlong <- msprep(time=c(NA,rep("time",6)),
            stat=c(NA,paste("stat",1:6,sep="")),
            data=ebmt2,keep=covs,trans=tmat)
# Actually run for a subset of the actual data
ebmtlongsub <- ebmtlong[1:1200,] # first 200 patients
# The reduced rank 2 solution
rr2 <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year+age, data=ebmtlongsub, R=2, eps=0.001)
rr2$Alpha; rr2$Gamma; rr2$Beta; rr2$loglik

# This reproduces the results in Fiocco, Putter & van Houwelingen (2005)
# Takes a while to run
\dontrun{
    # The reduced rank 3 solution
    rr3 <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year+age, data=ebmtlong, R=3)
    rr3$Alpha; rr3$Gamma; rr3$Beta; rr3$loglik
    # The reduced rank 3 solution, with no reduction on age
    rr3 <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year, full=~age, data=ebmtlong, R=3)
    rr3$Alpha; rr3$Gamma; rr3$Beta; rr3$loglik
    # The full rank solution
    fullrank <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year+age, data=ebmtlong, R=6)
    fullrank$Beta; fullrank$loglik
}
}
\keyword{survival}
