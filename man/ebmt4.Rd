\name{EBMT data}

\docType{data}

\alias{ebmt4}

\title{Data from the European Registry for Blood and Marrow Transplantation (EBMT)}

\description{A data frame of 2279 patients transplanted at the EBMT
between 1985 and 1998. These data were used in Fiocco, Putter & van Houwelingen (2008)
and van Houwelingen & Putter (2008).
The included variables are
  \describe{
  \item{id}{Patient identification number}
  \item{plat}{Time in days from transplantation to platelet recovery
  or last follow-up}
  \item{plat.s}{Platelet recovery status; 1 = platelet recovery, 0 =
  censored}
  \item{agvh}{Time in days from transplantation to acute graft-versus-host
  disease (AGvHD) or last follow-up}
  \item{agvh.s}{Acute GvHD status; 1 = acute GvHD, 0 = censored}
  \item{rel}{Time in days from transplantation to relapse or last follow-up}
  \item{rel.s}{Relapse status; 1 = relapse, 0 = censored}
  \item{srv}{Time in days from transplantation to death or last follow-up}
  \item{srv.s}{Relapse status; 1 = dead, 0 = censored}
  \item{year}{Year of transplantation; factor with levels "1985-1989",
  "1990-1994", "1995-1998"}
  \item{agecl}{Patient age at transplant; factor with levels "<=20",
  "20-40", ">40"}
  \item{tcd}{T-cell depletion; factor with levels "No TCD", "TCD"}
  \item{match}{Donor-recipient gender match; factor with levels
  "No gender mismatch", "Gender mismatch"}
  }
}

\usage{data(ebmt4)}

\format{A data frame, see \code{\link{data.frame}}.}

\references{
Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank proportional hazards regression
and simulation-based prediction for multi-state models. \emph{Statistics in Medicine} \bold{27},
4340--4358.

van Houwelingen HC, Putter H (2008). Dynamic predicting by landmarking as an alternative
for multi-state modeling: an application to acute lymphoid leukemia data. \emph{Lifetime Data Anal}
\bold{14}, 447--463.
}

\keyword{datasets}
