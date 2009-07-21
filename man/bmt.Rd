\name{bmt}
\alias{bmt}
\non_function{}
\title{BMT data from Klein and Moeschberger}
\description{
  A data frame of 137 rows (patients) and 22 columns. The included variables are
  \describe{
    \item{group}{
      Disease group; 1 = ALL, 2 = AML Low Risk, 3 = AML High Risk
    }
    \item{t1}{
      Time in days to death or last follow-up
    }
    \item{t2}{
      Disease-free survival time in days (time to relapse, death or last follow-up)
    }
    \item{d1}{
      Death indicator; 1 = dead, 0 = alive
    }
    \item{d2}{
      Relapse indicator; 1 = relapsed, 0 = disease-free
    }
    \item{d3}{
      Disease-free survival indicator; 1 = dead or relapsed, 0 = alive and disease-free)
    }
    \item{ta}{
      Time in days to Acute Graft-versus-Host Disease (AGVHD)
    }
    \item{da}{
      Acute GVHD indicator; 1 = Acute GVHD, 0 = No Acute GVHD
    }
    \item{tc}{
      Time (days) to Chronic Graft-vrsus-Host Disease (CGVHD)
    }
    \item{dc}{
      Chronic GVHD indicator; 1 = Chronic GVHD, 0 = No Chronic GVHD 
    }
    \item{tp}{
      Time (days) to platelet recovery
    }
    \item{dp}{
      Platelet recovery indicator; 1 = platelets returned to normal,
      0 = platelets never returned to normal
    }
    \item{z1}{
      Patient age in years
    }
    \item{z2}{
      Donor age in years
    }
    \item{z3}{
      Patient sex; 1 = male, 0 = female
    }
    \item{z4}{
      Donor sex; 1 = male, 0 = female
    }
    \item{z5}{
      Patient CMV status; 1 = CMV positive, 0 = CMV negative
    }
    \item{z6}{
      Donor CMV status; 1 = CMV positive, 0 = CMV negative
    }
    \item{z7}{
      Waiting time to transplant in days
    }
    \item{z8}{
      FAB; 1 = FAB grade 4 or 5 and AML, 0 = Otherwise
    }
    \item{z9}{
      Hospital; 1 = The Ohio State University, 2 = Alferd ,
      3 = St. Vincent, 4 = Hahnemann
    }
    \item{z10}{
      MTX used as a Graft-versus-Host prophylactic; 1 = yes, 0 = no
    }
  }
}
\format{
A data frame, see \code{\link{data.frame}}.
}
\references{
  Klein and Moeschberger (1997). \emph{Survival Analysis Techniques for Censored
    and Truncated Data}, Springer, New York.
}
\usage{
data(bmt)
}
\keyword{datasets}
