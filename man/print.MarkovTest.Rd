\name{print.MarkovTest}
\alias{print.MarkovTest}
\title{Print method for a MarkovTest object}
\description{Print method for an object of class 'MarkovTest'}
\usage{
\S3method{print}{MarkovTest}(x, ...)
}
\arguments{
  \item{x}{Object of class 'markovTest', as obtained by call to
  \code{\link{MarkovTest}}}
  \item{\dots}{Further arguments to print}
}
\value{No return value}
\seealso{\code{\link{MarkovTest}}}
\author{Hein Putter \email{H.Putter@lumc.nl}}
\examples{
# Example provided by the prothrombin data
data("prothr")
# Apply Markov test to grid of monthly time points over the first 7.5 years
year <- 365.25
month <- year / 12
grid <- month * (1 : 90)
# Markov test for transition 1 (wild bootstrap based on 100 replications)
MT <- MarkovTest(prothr, id = "id", transition = 1,
                 grid = grid, B = 100)
MT
}
\keyword{hplot}
