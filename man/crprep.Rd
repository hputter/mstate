\name{crprep}
\alias{crprep}
\title{Function to create weighted data set for competing risks analyses}
\description{
This function converts a dataset which is in short format (one subject per
line, one column indicating type of end point at end of follow-up) into a
weighted dataset in counting process style notation. With this data set,
competing risks analyses can be performed.
}
\usage{
crprep(Tstop, status, data, trans = 1, cens = 0, Tstart=0, id, keep,
       shorten = TRUE, rm.na = TRUE, origin = 0, prec.factor = 100)
}
\arguments{
  \item{Tstop}{Either 1) a vector containing the time at which the follow-up is
    ended, or 2) a character string indicating the column name that contains the
    end times. Missing values are allowed.}
  \item{status}{Either 1) a vector describing status at end of follow-up, or 2)
    a character string indicating the column name that contains this information.
    Missing values are allowed. See "Details".}
  \item{data}{Data frame in which to interpret "Tstart", "status", "trans",
    "cens", "Tstart", "id" and "keep", if appropriate.}
  \item{trans}{Values of "status" that are the event types of interest.
    Defaults to 1. See "Details".}
  \item{cens}{Value of "status" indicating censoring. Defaults to 0.}
  \item{Tstart}{Either 1) a vector containing the individual times at which the
    follow-up is started, or 2) a character string indicating the column name
    that contains the entry times, or 3) One numeric value in case it is the
    same for every subject. Missing values are allowed. Defaults to 0.}
  \item{id}{Either 1) a vector containing the subject identifiers, or 2) a
    character string indicating the column name containing these subject
    identifiers. If not provided, "id" will be assigned with values 1,...,n.}
  \item{keep}{Either 1) a data frame or matrix with n rows or a numeric or
    factor vector of length n containing covariate(s) that need to be retained
    in the output dataset, or 2) a character vector containing the column names
    of these covariates in data. See "Details".}
  \item{shorten}{Logical. If true, number of rows in output is reduced by
    collapsing rows within a subject in which weights do not change.}
  \item{rm.na}{Logical. If true, rows for which any of "Tstart", "status" or
    "Tstart" is missing are deleted.}
  \item{origin}{Substract origin time units from all Tstop and Tstart times.}
  \item{prec.factor}{Factor by which to multiply the machine's precision.
    Censoring and truncation times are shifted by prec.factor*precision if event
    times and censoring/trunaction times are equal.}
}
\details{
The function creates a data set that allows to perform analyses that are based
on the subdistribution hazard. For each type of outcome as specified via "trans",
individuals with a competing event remain in the risk set with a weight that
depends on the censoring and truncation mechanisms in the data. Typically, their
weights change over follow-up, and therefore such individuals are split into
several rows. Several types of outcome can be specified at once, thus allowing
for regression analyses using the "long format" data set (see Putter et al. 2007).
A regression on the cause-specific hazard using the created data set can be
performed by using "subset=count==0". 

If keep is a data.frame or a named matrix, the same names are used for the
covariate columns in the output data set. If keep is a matrix without names,
then the covariate columns are given the names "V1" until "Vk". If keep is a
vector from a (sub)-list, e.g. obj$name2$name1, then the column name is based
on the most inner part (i.e.\ "name1"). If keep is a vector of the form
obj[,"name1"], then the column is named "name1". For all other specifications,
the name is copied as is.

The current function does not allow to create a weighted data set in which the
censoring and /or truncation mechanisms depend on covariates. One option is to
create a weighted data set for each level of a categorical covariate.
}
\value{
A data frame in long (counting process) format containing the covariates
(replicated per subject) and the following columns 
  \item{patid}{subject identifier (1:n if argument "id" was missing)}
  \item{Tstart}{start dates of dataset (counting process notation)}
  \item{Tstop}{stop dates of dataset (counting process notation)}
  \item{status}{status of the subject at the end of his follow-up}
  \item{weight.cens}{weights due to censoring mechanism}
  \item{weight.trunc}{weights due to truncation mechanism (if present)}
  \item{count}{counter per subject and type of end point, 1 to number of rows
    per subject id and type of end point}
  \item{failcode}{type of end point, thus allowing to perform regression using 
    the "long format" data set}
  }
\references{
Geskus RB (2011). Cause-Specific Cumulative Incidence Estimation and the Fine
and Gray Model Under Both Left Truncation and Right Censoring.
\emph{Biometrics} \bold{67}, 39--49. 

Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine}
\bold{26}, 2389--2430.
}
\author{Ronald Geskus \email{statistics@inter.nl.net}, Kristian van Hemert;
  some adaptation by Marcel Wolbers}
\examples{
data(ebmt2)
Data <- ebmt2[c(1:5,7:8,30),] # just a small (non-representative) subset
Data$Tstart <- 0
Data$Tstop <- Data$time
Data$status[Data$status==6] <- 2 # recode 6 to 2
Data$dissub <- factor(Data$dissub) # get rid of missing level
Data.weight <- crprep(Data$Tstop, Data$status, trans=c(1,2), keep=Data$dissub)

# calculate cause-specific cumulative incidence, no truncation,
# compare with Cuminc (also from mstate)
ci <- Cuminc(Data$Tstop, Data$status)
ci
sf <- survfit(Surv(Tstart,Tstop,status==1)~1, data=Data.weight,
              weight=weight.cens, subset=failcode==1)
sf <- summary(sf)
data.frame(time=sf$time, CI.1=1-sf$surv)
sf <- survfit(Surv(Tstart,Tstop,status==2)~1, data=Data.weight,
              weight=weight.cens, subset=failcode==2)
sf <- summary(sf)
data.frame(time=sf$time, CI.2=1-sf$surv)

# Fine and Gray regression for cause 1
cw <- coxph(Surv(Tstart,Tstop,status==1)~dissub, data=Data.weight,
      weight=weight.cens, subset=failcode==1)
cw
# This can be checked with the results of crr (cmprsk)
# crr(ftime=Data$Tstop, fstatus=Data$status, cov1=as.numeric(Data$dissub)-1)

\donttest{
# Create simulated data
set.seed(1234)
N <- 200
p <- 0.3
p.t <- 0.6
Z1 <- rnorm(N)
Z2 <- rnorm(N)
tmp <- runif(N)
Data <- data.frame(Tstart=runif(N,0,3)*rbinom(N,1,p.t),
          Tstop=ifelse( tmp>(1-p)^(exp(0.5*Z1+0.5*Z2)),
                  -log(1-(1-exp(log(tmp)/exp(0.5*Z1+0.5*Z2)))/p),
                  -log(tmp)/exp(-0.5*Z1+0.5*Z2)),
          stat=ifelse(tmp>(1-p)^(exp(0.5*Z1+0.5*Z2)),1,2),
          Z1=Z1, Z2=Z2, tstat=rep(0,N))
Data$type <- Data$stat
Data$tevent <- Data$Tstop
Data$cens <- runif(N,0.5,1)
Data[Data$cens<Data$Tstop,"stat"] <- 0
Data$Tstop <- pmin(Data$Tstop,Data$cens)
Data$tstat <- ifelse(Data$Tstart > Data$Tstop, 1, 0)
Data.weight <- crprep(Data$Tstop, Data$stat, trans=c(1,2))

# calculate cause-specific cumulative incidence, no truncation,
# compare with Cuminc (also from mstate)
ci <- Cuminc(Data$Tstop, Data$stat)
plot(ci$time,ci$CI.1,type="s",lwd=3,col="black",ylim=c(0,0.5),
   xlab="Time",ylab="Cumulative incidence")
lines(ci$time,ci$CI.2,type="s",lwd=3,col="red")
lines(ci$time,ci$CI.1-qnorm(0.975)*ci$seCI.1,type="s",lty=3)
lines(ci$time,ci$CI.1+qnorm(0.975)*ci$seCI.1,type="s",lty=3)
lines(survfit(Surv(Tstart,Tstop,status==1)~1,data=Data.weight,
   weight=weight.cens,subset=failcode==1),
   fun="event",col="lightblue",lwd=1,mark.time=FALSE)
lines(survfit(Surv(Tstart,Tstop,status==1)~1,data=Data.weight,
   weight=weight.cens,subset=failcode==1),
   fun="event",col="lightblue",mark.time=FALSE,conf.int="only",lty=2)


# Proportional hazards regression on subdistribution and cause-specific hazard,
# with truncation
Data.weight.trunc <- crprep("Tstop", "stat", Tstart="Tstart",
                       data=subset(Data,tstat==0), trans=1:2, keep=c("Z1","Z2"))
coxph(Surv(Tstart,Tstop,status==1)~Z1+Z2,data=Data.weight.trunc,
   weight=weight.cens*weight.trunc,subset=failcode==1) #cause 1
# Both end points, assume effect Z2 same for both
coxph(Surv(Tstart,Tstop,status==1)~strata(failcode)*Z1+Z2,
   data=Data.weight.trunc,weight=weight.cens*weight.trunc)
# Cause-specific hazard
coxph(Surv(Tstart,Tstop,stat==1)~Z1+Z2,data=subset(Data, tstat==0))
coxph(Surv(Tstart,Tstop,status==1)~Z1+Z2,data=Data.weight.trunc,
   subset=failcode==1&count==1)

data(ebmt2)
ebmt2.long <- crprep("time", "status", data=ebmt2, trans=1:6,
                keep=c("dissub", "match", "tcd", "year", "age"))
# ebmt2.long <- with(ebmt2, crprep(time, cod, trans=levels(cod)[-1], cens="Alive",
#                  keep=c("dissub", "match", "tcd", "year", "age")))
plot(survfit(Surv(time, status>0)~1, data=ebmt2, etype=cod), mark.time=FALSE,
   col=2:7, fun="event",lwd=3)
lines(survfit(Surv(Tstart, Tstop, failcode==status) ~ failcode, data=ebmt2.long,
   weight=weight.cens), lwd=2, lty=2, fun="event", mark.time=FALSE)
}
}
\keyword{datagen}
