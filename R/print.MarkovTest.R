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
