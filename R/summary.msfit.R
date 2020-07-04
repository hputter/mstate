summary.msfit <- function(object, times, transitions,
                          variance=TRUE, conf.int=0.95,
                          conf.type=c("log", "none", "plain"),
                          extend=FALSE, ...)
{
    if (!inherits(object, "msfit"))
        stop("'object' must be a 'msfit' object")

    conf.type <- match.arg(conf.type)
    if (!conf.type %in% c("log", "none", "plain"))
        stop("conf.type should be one of log, none, plain")
    Haz <- object$Haz
    varHaz <- object$varHaz
    K <- max(Haz$trans)
    tt <- unique(Haz$time) # the time points
    nt <- length(tt)
    
    if (missing(transitions)) transitions <- 1:K
    transitions <- intersect(transitions, 1:K)
    
    res <- list()
    if (variance) {
        if (missing(times)) {
            for (k in transitions) {
                dfr <- Haz[Haz$trans == k, ]
                varHazk <- varHaz$varHaz[varHaz$trans1 == k & varHaz$trans2 == k]
                dfr[, 3] <- sqrt(varHazk)
                names(dfr)[3] <- "seHaz"
                # Calculate confidence intervals
                if (!conf.type=="none") {
                    if (conf.type=="plain") {
                        lower <- dfr$Haz - qnorm(1 - (1-conf.int)/2) * dfr$seHaz
                        upper <- dfr$Haz + qnorm(1 - (1-conf.int)/2) * dfr$seHaz
                        lower[lower<0] <- 0
                    }
                    else if (conf.type=="log") {
                        lower <- exp(log(dfr$Haz) - qnorm(1 - (1-conf.int)/2) * dfr$seHaz / dfr$Haz)
                        upper <- exp(log(dfr$Haz) + qnorm(1 - (1-conf.int)/2) * dfr$seHaz / dfr$Haz)
                        lower[lower<0] <- 0
                        lower[dfr$seHaz==0] <- dfr$Haz[dfr$seHaz==0]
                        upper[dfr$seHaz==0] <- dfr$Haz[dfr$seHaz==0]
                    }
                    dfr <- cbind(dfr, lower, upper)
                    names(dfr)[4:5] <- c("lower", "upper")
                }
                res[[k]] <- dfr
            }
        } else {
            if (!extend) times <- times[times<=max(tt)]
            for (k in transitions) {
                Hazk <- Haz[Haz$trans == k, ]
                varHazk <- varHaz$varHaz[varHaz$trans1 == k & varHaz$trans2 == k]
                dfr <- matrix(NA, length(times), 3)
                dfr[, 1] <- times
                dfr[, 2] <- approx(x=tt, y=Hazk$Haz, xout=times,
                                   f=0, method="constant", rule=2)$y
                dfr[, 3] <- approx(x=tt, y=sqrt(varHazk), xout=times,
                                   f=0, method="constant", rule=2)$y
                dfr <- as.data.frame(dfr)
                names(dfr) <- c("times", "Haz", "seHaz")
                # Calculate confidence intervals
                if (!conf.type=="none") {
                    if (conf.type=="plain") {
                        lower <- dfr$Haz - qnorm(1 - (1-conf.int)/2) * dfr$seHaz
                        upper <- dfr$Haz + qnorm(1 - (1-conf.int)/2) * dfr$seHaz
                        lower[lower<0] <- 0
                    }
                    else if (conf.type=="log") {
                        lower <- exp(log(dfr$Haz) - qnorm(1 - (1-conf.int)/2) * dfr$seHaz / dfr$Haz)
                        upper <- exp(log(dfr$Haz) + qnorm(1 - (1-conf.int)/2) * dfr$seHaz / dfr$Haz)
                        lower[lower<0] <- 0
                        lower[dfr$seHaz==0] <- dfr$Haz[dfr$seHaz==0]
                        upper[dfr$seHaz==0] <- dfr$Haz[dfr$seHaz==0]
                    }
                    dfr <- cbind(dfr, lower, upper)
                    names(dfr)[4:5] <- c("lower", "upper")
                }
                res[[k]] <- dfr
            }
        }
    } else { # no variance
        if (missing(times)) {
            for (k in transitions)
                res[[k]] <- Haz[Haz$trans == k, c("time", "Haz")]
        } else {
            if (!extend) times <- times[times<=max(tt)]
            for (k in transitions) {
                Hazk <- Haz[Haz$trans == k, ]
                dfr <- matrix(NA, length(times), 2)
                dfr[, 1] <- times
                dfr[, 2] <- approx(x=tt, y=Hazk$Haz, xout=times,
                                   f=0, method="constant", rule=2)$y
                dfr <- as.data.frame(dfr)
                names(dfr) <- c("times", "Haz")
                res[[k]] <- dfr
            }
        }
    }
    attr(res, "transitions") <- transitions
    class(res) <- "summary.msfit"
    return(res)
}

print.summary.msfit <- function(x, complete=FALSE, ...)
{
    if (!inherits(x, "summary.msfit"))
        stop("'x' must be a 'summary.msfit' object")
    transitions <- attr(x, "transitions")
    tt <- unique(x[[1]]$time) # the time points
    nt <- length(tt)
    if (nt<=12 | complete) {
        for (k in transitions) {
            cat("\nTransition", k, ":\n")
            ptk <- x[[k]]
            print(ptk, ...)
        }
    } else {
        for (k in transitions) {
            cat("\nTransition", k, "(head and tail):\n")
            ptk <- x[[k]]
            print(head(ptk), ...)
            cat("\n...\n")
            print(tail(ptk), ...)
        }
    }
    return(invisible())
}
