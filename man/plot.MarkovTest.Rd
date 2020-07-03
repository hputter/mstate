\name{plot.MarkovTest}
\alias{plot.MarkovTest}
\title{Plot method for a MarkovTest object}
\description{
Plot method for an object of class 'MarkovTest'. It plots the trace of the log-rank
statistics provided by \code{\link{MarkovTest}}.}
\usage{
\S3method{plot}{MarkovTest}(x, y,...)
}
\arguments{
  \item{x}{Object of class 'probtrans', containing estimated transition
  probabilities}
  \item{\dots}{Further arguments to plot}
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

plot(MT, grid, what="states", idx=1:50, states=rownames(attr(prothr, "trans")),
     xlab="Days since randomisation", ylab="Log-rank test statistic",
     main="Transition Normal -> Low")

plot(MT, grid,what="overall", idx=1:50,
     xlab="Days since randomisation", ylab="Chi-square test statistic",
     main="Transition Normal -> Low")

plot(MT, grid, what="states", quantiles=FALSE) # only trace
plot(MT, grid, what="states") # trace plus quantiles (default)
plot(MT, grid, what="states", idx=1:10) # trace plus quantiles, plus first 10 bootstrap traces
 
plot(MT, grid, what="overall", quantiles=FALSE) # only trace
plot(MT, grid, what="overall") # trace plus quantiles (default)
plot(MT, grid, what="overall", idx=1:10) # trace plus quantiles, plus first 10 bootstrap traces
}
\keyword{hplot}
