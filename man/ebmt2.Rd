\name{EBMT cause of death data}
\docType{data}
\alias{ebmt2}
\title{Data from the European Society for Blood and Marrow Transplantation (EBMT)}
\description{A data frame of 8966 patients transplanted at the EBMT.
The included variables are
  \describe{
  \item{id}{Patient identification number}
  \item{time}{Time in months from transplantation to death or last
  follow-up}
  \item{status}{Survival status; 0 = censored; 1,...,6 = death
  due to the following causes: Relapse (1), GvHD (2),
  Bacterial infections (3), Viral infections (4),
  Fungal infections (5), Other causes (6)}
  \item{cod}{Cause of death as factor with levels "Alive", "Relapse",
  "GvHD", "Bacterial", "Viral", "Fungal", "Other"}
  \item{dissub}{Disease subclassification; factor with levels "AML",
  "ALL", "CML"}
  \item{match}{Donor-recipient gender match; factor with levels
  "No gender mismatch", "Gender mismatch"}
  \item{tcd}{T-cell depletion; factor with levels "No TCD", "TCD",
  "Unknown"}
  \item{year}{Year of transplantation; factor with levels
  "1985-1989", "1990-1994", "1995-1998"}
  \item{age}{Patient age at transplant; factor with levels "<=20",
  "20-40", ">40"}
  }
}
\usage{data(ebmt2)}
\format{A data frame, see \code{\link{data.frame}}.}
\source{We acknowledge the European Society for Blood and Marrow Transplantation
(EBMT) for making available these data. Disclaimer: these data were simplified
for the purpose of illustration of the analysis of competing risks and multi-state
models and do not reflect any real life situation. No clinical conclusions should
be drawn from these data.}
\references{
Fiocco M, Putter H, van Houwelingen JC (2005). Reduced rank proportional
hazards model for competing risks. \emph{Biostatistics} \bold{6}, 465--478.
}
\keyword{datasets}
