#' Summary method for a probtrans object
#' 
#' Summary method for an object of class 'probtrans'. It prints a selection of
#' the estimated transition probabilities, and, if requested, also of the
#' variances.
#' 
#' 
#' @aliases summary.probtrans print.summary.probtrans
#' @param object Object of class 'probtrans', containing estimated transition
#' probabilities from and to all states in a multi-state model
#' @param times Time points at which to evaluate the transition probabilites
#' @param from Specifies from which state the transition probabilities are to
#' be printed. Should be subset of 1:S, with S the number of states in the
#' multi-state model. Default is print from state 1 only. User can specify
#' from=0 to print transition probabilities from all states
#' @param to Specifies the transition probabilities to which state are to be
#' printed. User can specify to=0 to print transition probabilities to all
#' states. This is also the default
#' @param variance Whether or not the standard errors of the estimated
#' transition probabilities should be printed; default is \code{TRUE}
#' @param conf.int The proportion to be covered by the confidence intervals,
#' default is 0.95
#' @param conf.type The type of confidence interval, one of "log", "none", or
#' "plain". Defaults to "log"
#' @param extend logical value: if \code{TRUE}, prints information for all
#' specified times, even if there are no subjects left at the end of the
#' specified times. This is only valid if the times argument is present
#' @param x Object of class 'summary.probtrans', to be printed
#' @param complete Whether or not the complete estimated transition
#' probabilities should be printed (\code{TRUE}) or not (\code{FALSE}); default
#' is \code{FALSE}, in which case the estimated transition probilities will be
#' printed for the first and last 6 time points of each starting state or of
#' the selected times (or all when there are at most 12 of these time points
#' @param \dots Further arguments to print
#' @return Function \code{summary.probtrans} returns an object of class
#' "summary.probtrans", which is a list (for each \code{from} state) of
#' transition probabilities at the specified (or all) time points. The
#' \code{print} method of a \code{summary.probtrans} doesn't return a value.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @seealso \code{\link{probtrans}}
#' @keywords print
#' @examples
#' 
#' # First run the example of probtrans
#' tmat <- trans.illdeath()
#' tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'                  dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'                  x1=c(1,1,1,0,0,0),x2=c(6:1))
#' tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
#'                  data=tg,keep=c("x1","x2"),trans=tmat)
#' tglong <- expand.covs(tglong,c("x1","x2"))
#' cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
#'             data=tglong,method="breslow")
#' newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
#' HvH <- msfit(cx,newdata,trans=tmat)
#' pt <- probtrans(HvH,predt=0)
#' 
#' # Default, prediction from state 1
#' summary(pt)
#' # Only from states 1 and 3
#' summary(pt, from=c(1, 3))
#' # Use from=0 for prediction from all states
#' summary(pt, from=0)
#' # Only to states 1 and 2
#' summary(pt, to=1:2)
#' # Default is 95% confidence interval, change here to 90% 
#' summary(pt, to=1:2, conf.int=0.90)
#' # Do not show variances (nor confidence intervals)
#' summary(pt, to=1:2, variance=FALSE)
#' # Transition probabilities only at specified time points
#' summary(pt, times=seq(0, 15, by=3))
#' # Last specified time point is larger than last observed, not printed
#' # Use extend=TRUE as in summary.survfit
#' summary(pt, times=seq(0, 15, by=3), extend=TRUE)
#' # Different types of confidence intervals, default is log
#' summary(pt, times=seq(0, 15, by=3), conf.type="plain")
#' summary(pt, times=seq(0, 15, by=3), conf.type="no")
#' # When the number of time points specified is larger than 12, head and tail is shown
#' x <- summary(pt, times=seq(5, 8, by=0.25))
#' x
#' # If all time points should be printed, specify complete=TRUE in the print statement
#' print(x, complete=TRUE)
#' 
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
