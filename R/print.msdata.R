#' Print method for a msdata object
#' 
#' Print method for an object of class 'msdata'
#' 
#' 
#' @param x Object of class 'msdata', as prepared for instance by
#' \code{\link{msprep}}
#' @param trans Boolean specifying whether or not the transition matrix should
#' be printed as well; default is \code{FALSE}
#' @param \dots Further arguments to print
#' @return No return value
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @seealso \code{\link{probtrans}}
#' @keywords hplot
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # some data in wide format
#' tg <- data.frame(stt=rep(0,6),sts=rep(0,6),
#'         illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'         dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'         x1=c(1,1,1,2,2,2),x2=c(6:1))
#' tg$x1 <- factor(tg$x1,labels=c("male","female"))
#' tg$patid <- factor(2:7,levels=1:8,labels=as.character(1:8))
#' # define time, status and covariates also as matrices
#' tt <- matrix(c(rep(NA,6),tg$illt,tg$dt),6,3)
#' st <- matrix(c(rep(NA,6),tg$ills,tg$ds),6,3)
#' keepmat <- data.frame(gender=tg$x1,age=tg$x2)
#' # data in long format using msprep
#' msp <- msprep(time=tt,status=st,trans=tmat,keep=as.matrix(keepmat))
#' print(msp)
#' print(msp, trans=TRUE)
#' 
print.msdata <- function(x,trans=FALSE,...)
{
    if (!inherits(x, "msdata"))
        stop("'x' must be an 'msdata' object")
    cat("An object of class 'msdata'\n\nData:\n")
    print.data.frame(x)
    if (trans) {
        trans <- attr(x, "trans")
        cat("\nTransition matrix:\n")
        print(trans)
    }
    return(invisible())
}
