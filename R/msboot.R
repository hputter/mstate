#' Bootstrap function in multi-state models
#' 
#' A generic nonparametric bootstrapping function for multi-state models.
#' 
#' The function \code{msboot} samples randomly with replacement subjects from
#' the original dataset \code{data}. The individuals are identified with
#' \code{id}, and bootstrap datasets are produced by concatenating all selected
#' rows.
#' 
#' @param theta A function of \code{data} and perhaps other arguments,
#' returning the value of the statistic to be bootstrapped; the output of theta
#' should be a scalar or numeric vector
#' @param data An object of class 'msdata', such as output from
#' \code{\link{msprep}}
#' @param B The number of bootstrap replications; the default is taken to be
#' quite small (5) since bootstrapping can be time-consuming
#' @param id Character string indicating which column identifies the subjects
#' to be resampled
#' @param verbose The level of output; default 0 = no output, 1 = print the
#' replication
#' @param ... Any further arguments to the function \code{theta}
#' @return Matrix of dimension (length of output of theta) x B, with b'th
#' column being the value of theta for the b'th bootstrap dataset
#' @author Marta Fiocco, Hein Putter <H.Putter@@lumc.nl>
#' @references Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
#' proportional hazards regression and simulation-based prediction for
#' multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.
#' @keywords datagen
#' @examples
#' 
#' tmat <- trans.illdeath()
#' data(ebmt1)
#' covs <- c("score","yrel")
#' msebmt <- msprep(time=c(NA,"rel","srv"),status=c(NA,"relstat","srvstat"),
#' 		data=ebmt1,id="patid",keep=covs,trans=tmat)
#' # define a function (this one returns vector of regression coef's)
#' regcoefvec <- function(data) {
#'   cx <- coxph(Surv(Tstart,Tstop,status)~score+strata(trans),
#'           data=data,method="breslow")
#'   return(coef(cx))
#' }
#' regcoefvec(msebmt)
#' set.seed(1234)
#' msboot(theta=regcoefvec,data=msebmt,id="patid")
#' 
#' @export msboot
`msboot` <- function(theta,data,B=5,id="id",verbose=0,...)
{
    if (!inherits(data, "msdata"))
        stop("'data' must be a 'msdata' object")
    trans <- attr(data, "trans")
    ids <- unique(data[[id]])
    n <- length(ids)
    th <- theta(data,...) # actually only used to get the length
    res <- matrix(NA,length(th),B)
    for (b in 1:B) {
        if (verbose>0) {
            cat("\nBootstrap replication",b,"\n")
            flush.console()
        }
        bootdata <- NULL
        bids <- sample(ids,replace=TRUE)
        bidxs <- unlist(sapply(bids, function(x) which(x==data[[id]])))
        bootdata <- data[bidxs,]
        if (verbose>0) {
            print(date())
            print(events(bootdata))
            cat("applying theta ...")
        }
        thstar <- theta(bootdata,...)
        res[,b] <- thstar
    }
    if (verbose) cat("\n")
    return(res)
}
