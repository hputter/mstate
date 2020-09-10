#' Print method for a MarkovTest object
#' 
#' Print method for an object of class 'MarkovTest'
#' 
#' 
#' @param x Object of class 'markovTest', as obtained by call to
#' \code{\link{MarkovTest}}
#' @param \dots Further arguments to print
#' @return No return value
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @seealso \code{\link{MarkovTest}}
#' @keywords hplot
#' @examples
#' 
#' # Example provided by the prothrombin data
#' data("prothr")
#' # Apply Markov test to grid of monthly time points over the first 7.5 years
#' year <- 365.25
#' month <- year / 12
#' grid <- month * (1 : 90)
#' # Markov test for transition 1 (wild bootstrap based on 25 replications for brevity)
#' MT <- MarkovTest(prothr, id = "id", transition = 1,
#'                  grid = grid, B = 25)
#' MT
#' 
#' @export
print.MarkovTest <- function(x, ...)
{
  cat("Log-rank based Markov test for transition ", x$trans, " (", x$from, " -> ", x$to, ")\n", sep="")
  cat("Qualifying states: ", x$qualset, "\n")
  cat("\nP-values based on ", x$B, " wild bootstrap replications with ")
  if (x$dist == "poisson") cat("centred Poisson")
  if (x$dist == "normal") cat("standard normal")
  cat(" distributed random weights\n")
  
  cat("\nQualifying state specific tests:\n")
  nfuns <- nrow(x$orig_stat)
  for (j in 1:nfuns) {
    cat("\n")
    print(x$fn[[1]][[j]])
    cat("\n")
    tmp <- t(rbind(x$qualset, x$orig_stat[j, ], x$p_stat_wb[j, ]))
    tmp <- as.data.frame(tmp)
    names(tmp) <- c("Qualstate", "u", "P(|U| > u)")
    print(tmp, row.names = FALSE)
  }
  
  cat("\nOverall test:\n")
  nfuns <- length(x$orig_ch_stat)
  for (j in 1:nfuns) {
    cat("\n")
    print(x$fn2[[j]])
    cat("\n")
    tmp <- matrix(c(x$orig_ch_stat[j], x$p_ch_stat_wb), 1, 2)
    tmp <- as.data.frame(tmp)
    names(tmp) <- c("chisq", "P(|U| > u)")
    print(tmp, row.names = FALSE)
  }
}
