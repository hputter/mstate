summary.probtrans <- function(object, times, from=1, to=0,
                              variance=TRUE, conf.int=0.95,
                              conf.type=c("log", "none", "plain"),
                              extend=FALSE, ...)
{
    if (!inherits(object, "probtrans"))
        stop("'object' must be a 'probtrans' object")
    conf.type <- match.arg(conf.type)
    if (!conf.type %in% c("log", "none", "plain"))
        stop("conf.type should be one of log, none, plain")
    trans <- object$trans
    S <- dim(trans)[1]
    tt <- unique(object[[1]]$time) # the time points
    nt <- length(tt)
    if (!all(from %in% 0:S)) stop("from should be either 0 (all states) or 1 to S (number of states)")
    if (!all(to %in% 0:S)) stop("to should be either 0 (all states) or 1 to S (number of states)")
    if (any(from == 0)) from <- 1:S
    if (any(to==0)) to <- 1:S
    cols <- 1 + to
    probcols <- cols
    if (variance & ncol(object[[1]]) == S+1) {
        warning("probtrans object does not contain standard errors, setting option variance to FALSE")
        variance <- FALSE
    }
    if (variance) {
        secols <- S + 1 + to
        cols <- c(cols, secols)
        if (missing(times)) {
            res <- list()
            for (k in from) {
                ptk <- object[[k]]
                res[[k]] <- ptk # Default in case of no confidence intervals
                if (!conf.type=="none") {
                    if (conf.type=="plain") {
                        lower <- ptk[, probcols] - qnorm(1 - (1-conf.int)/2) * ptk[, secols]
                        upper <- ptk[, probcols] + qnorm(1 - (1-conf.int)/2) * ptk[, secols]
                        lower[lower<0] <- 0
                        upper[upper>1] <- 1
                    }
                    else if (conf.type=="log") {
                        lower <- exp(log(ptk[, probcols]) - qnorm(1 - (1-conf.int)/2) * ptk[, secols] / ptk[, probcols])
                        upper <- exp(log(ptk[, probcols]) + qnorm(1 - (1-conf.int)/2) * ptk[, secols] / ptk[, probcols])
                        lower[lower<0] <- 0
                        upper[upper>1] <- 1
                        lower[ptk[, secols]==0] <- ptk[, probcols][ptk[, secols]==0]
                        upper[ptk[, secols]==0] <- ptk[, probcols][ptk[, secols]==0]
                    }
                    dfr <- cbind(ptk$time, ptk[, c(probcols, secols)], lower, upper)
                    names(dfr)[1] <- "time"
                    names(dfr)[2 * length(to) + 1 + (1:length(to))] <- paste("lower", to, sep="")
                    names(dfr)[3 * length(to) + 1 + (1:length(to))] <- paste("upper", to, sep="")
                } else {
                    dfr <- cbind(ptk$time, ptk[, c(probcols, secols)])
                    names(dfr)[1] <- "time"
                }
                res[[k]] <- dfr
            }
        } else {
            if (!extend) {
                times <- times[times<=max(tt)]
            }
            res <- list()
            for (k in from) {
                ptk <- object[[k]]
                dfr <- matrix(NA, length(times), 1 + length(cols))
                dfr[, 1] <- times
                for (i in 1:length(cols))
                    dfr[, 1 + i] <- approx(x=tt, y=ptk[, cols[i]], xout=times,
                                           f=0, method="constant", rule=2)$y
                dfr <- as.data.frame(dfr)
                names(dfr) <- c("times", names(ptk)[cols])
                probcols <- 1 + (1:length(to))
                secols <- length(to) + 1 + (1:length(to))
                # confidence intervals uitrekenen
                if (!conf.type=="none") {
                    if (conf.type=="plain") {
                        lower <- dfr[, probcols] - qnorm(1 - (1-conf.int)/2) * dfr[, secols]
                        upper <- dfr[, probcols] + qnorm(1 - (1-conf.int)/2) * dfr[, secols]
                        lower[lower<0] <- 0
                        upper[upper>1] <- 1
                    }
                    else if (conf.type=="log") {
                        lower <- exp(log(dfr[, probcols]) - qnorm(1 - (1-conf.int)/2) * dfr[, secols] / dfr[, probcols])
                        upper <- exp(log(dfr[, probcols]) + qnorm(1 - (1-conf.int)/2) * dfr[, secols] / dfr[, probcols])
                        lower[lower<0] <- 0
                        upper[upper>1] <- 1
                        lower[dfr[, secols]==0] <- dfr[, probcols][dfr[, secols]==0]
                        upper[dfr[, secols]==0] <- dfr[, probcols][dfr[, secols]==0]
                    }
                    dfr <- cbind(dfr, lower, upper)
                    names(dfr)[2 * length(to) + 1 + (1:length(to))] <- paste("lower", to, sep="")
                    names(dfr)[3 * length(to) + 1 + (1:length(to))] <- paste("upper", to, sep="")
                }
                res[[k]] <- dfr
            }
        }
    }
    else {
        if (missing(times)) {
            res <- list()
            for (k in from) {
                ptk <- object[[k]]
                res[[k]] <- ptk # Default in case of no confidence intervals
                dfr <- cbind(ptk$time, ptk[, probcols])
                names(dfr)[1] <- "time"
                res[[k]] <- dfr
            }
        } else {
            if (!extend) {
                times <- times[times<=max(tt)]
            }
            res <- list()
            for (k in from) {
                ptk <- object[[k]]
                dfr <- matrix(NA, length(times), 1 + length(cols))
                dfr[, 1] <- times
                for (i in 1:length(cols))
                    dfr[, 1 + i] <- approx(x=tt, y=ptk[, cols[i]], xout=times,
                                           f=0, method="constant", rule=2)$y
                dfr <- as.data.frame(dfr)
                names(dfr) <- c("times", names(ptk)[cols])
                probcols <- 1 + (1:length(to))
                res[[k]] <- dfr
            }
        }
    }
    attr(res, "from") <- from
    class(res) <- "summary.probtrans"
    return(res)
}

print.summary.probtrans <- function(x, complete=FALSE, ...)
{
    if (!inherits(x, "summary.probtrans"))
        stop("'x' must be a 'summary.probtrans' object")
    from <- attr(x, "from")
    tt <- unique(x[[1]]$time) # the time points
    nt <- length(tt)
    if (nt<=12 | complete) {
        for (k in from) {
            cat("\nPrediction from state", k, ":\n")
            ptk <- x[[k]]
            print(ptk, ...)
        }
    } else {
        for (k in from) {
            cat("\nPrediction from state", k, "(head and tail):\n")
            ptk <- x[[k]]
            print(head(ptk), ...)
            cat("\n...\n")
            print(tail(ptk), ...)
        }
    }
    return(invisible())
}
