#' Summary method for an msfit object
#' 
#' Summary method for an object of class 'msfit'. It prints a selection of the
#' estimated cumulative transition intensities, and, if requested, also of the
#' (co)variances.
#' 
#' 
#' @aliases summary.msfit print.summary.msfit
#' @param object Object of class 'msfit', containing estimated cumulative
#' transition intensities for all transitions in a multi-state model
#' @param times Time points at which to evaluate the cumulative transition
#' hazards
#' @param transitions The transition for which to summarize the cumulative
#' transition hazards
#' @param variance Whether or not the standard errors of the estimated
#' cumulative transition intensities should be printed; default is \code{TRUE}
#' @param conf.int The proportion to be covered by the confidence intervals,
#' default is 0.95
#' @param conf.type The type of confidence interval, one of "log", "none", or
#' "plain". Defaults to "log"
#' @param extend logical value: if \code{TRUE}, prints information for all
#' specified times, even if there are no subjects left at the end of the
#' specified times. This is only valid if the times argument is present
#' @param x Object of class 'summary.probtrans', to be printed
#' @param complete Whether or not the complete estimated cumulative transition
#' intensities should be printed (\code{TRUE}) or not (\code{FALSE}); default
#' is \code{FALSE}, in which case the estimated cumulative transition hazards
#' will be printed for the first and last 6 time points of each transition or
#' of the selected times (or all when there are at most 12 of these time points
#' @param \dots Further arguments to summary
#' @return Function \code{summary.msfit} returns an object of class
#' "summary.msfit", which is a list (for each \code{from} state) of cumulative
#' transition hazaards at the specified (or all) time points. The \code{print}
#' method of a \code{summary.probtrans} doesn't return a value.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @seealso \code{\link{msfit}}
#' @keywords print
#' @examples
#' 
#' # Start with example from msfit
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
#' msf <- msfit(cx,newdata,trans=tmat)
#' 
#' # Default, all transitions, with SE
#' summary(msf)
#' summary(msf, conf.type="plain")
#' # Only transitions 1 and 3
#' summary(msf, tra=c(1, 3))
#' # Default is 95% confidence interval, change here to 90% 
#' summary(msf, conf.int=0.90)
#' # Do not show variances (nor confidence intervals)
#' summary(msf, variance=FALSE)
#' # Cumulative hazards only at specified time points
#' summary(msf, times=seq(0, 15, by=3))
#' # Last specified time point is larger than last observed, not printed
#' # Use extend=TRUE as in summary.survfit
#' summary(msf, times=seq(0, 15, by=3), extend=TRUE)
#' # Different types of confidence intervals, default is log
#' summary(msf, times=seq(0, 15, by=3), conf.type="plain")
#' summary(msf, times=seq(0, 15, by=3), conf.type="no")
#' # When the number of time points specified is larger than 12, head and tail is shown
#' x <- summary(msf, times=seq(5, 8, by=0.25))
#' x
#' # If all time points should be printed, specify complete=TRUE in the print statement
#' print(x, complete=TRUE)
#' 
#' @method summary msfit
#' @export
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
