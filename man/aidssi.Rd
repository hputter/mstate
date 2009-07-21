\name{HIV/AIDS SI data}

\docType{data}

\alias{aidssi}

\title{Data from the Amsterdam Cohort Studies on HIV infection and AIDS}

\description{A data frame of 329 homosexual men from the Amsterdam Cohort Studies
on HIV infection and AIDS. The included variables are
  \describe{
  \item{patid}{Patient identification number}
  \item{time}{Time in years from HIV infection to either SI appearance, AIDS, or
  last follow-up}
  \item{status}{Event indicator; 0 = censored, 1 = AIDS, 2 = SI appearance}
  \item{cause}{Failure cause; factor with levels "event-free", "AIDS", "SI"}
  \item{ccr5}{CCR5 genotype; factor with levels "WW" (2 wild-type alleles),
  "WM" (1 (or possibly two) mutant alleles)}
  }
}

\usage{data(aidssi)}

\format{A data frame, see \code{\link{data.frame}}.}

\details{This data was used in the tutorial on competing risks and multi-state
models (Putter, Fiocco, Geskus, 2007). For more information refer to the tutorial
or to the original papers (Geskus et al., 2000, 2003).}

\references{
Geskus RB (2000). On the inclusion of prevalent cases in HIV/AIDS natural
history studies through a marker-based estimate of time since seroconversion.
\emph{Statistics in Medicine} \bold{19}, 1753--1769.

Geskus RB, Miedema FA, Goudsmit J, Reiss P, Schuitemaker H, Coutinho RA (2003).
Prediction of residual time to AIDS and death based on markers and cofactors.
\emph{Journal of AIDS} \bold{32}, 514--521.

Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
2389--2430.
}

\keyword{datasets}
