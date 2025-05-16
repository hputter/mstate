#' Summary method for a probtrans.subjects object
#' 
#' Summary method for an object of class 'probtrans.subjects'. It prints a selection of
#' the estimated transition probabilities. Wrapper for 
#' \code{\link{summary.probtrans}}.
#' 
#' @aliases summary.probtrans.subjects
#' @param object Object of class 'probtrans.subjects', containing estimated transition
#' probabilities from and to all states in a multi-state model
#' @param id Subject identifier
#' @param times Time points at which to evaluate the transition probabilites
#' @param from Specifies from which state the transition probabilities are to
#' be printed. Should be subset of 1:S, with S the number of states in the
#' multi-state model. Default is print from state 1 only. User can specify
#' from=0 to print transition probabilities from all states
#' @param to Specifies the transition probabilities to which state are to be
#' printed. User can specify to=0 to print transition probabilities to all
#' states. This is also the default
#' @param extend logical value: if \code{TRUE}, prints information for all
#' specified times, even if there are no subjects left at the end of the
#' specified times. This is only valid if the times argument is present
#' @param \dots Further arguments to \code{\link{summary.probtrans}}
#' 
#' 
#' @import mstate
#' 
#' @return Function \code{summary.probtrans} returns an object of class
#' "summary.probtrans", which is a list (for each \code{from} state) of
#' transition probabilities at the specified (or all) time points. The
#' \code{print} method of a \code{summary.probtrans} doesn't return a value.
#' @author Hein Putter and Daniel Gomon
#' @seealso \code{\link[mstate:summary.probtrans]{summary.probtrans}}
#' @keywords print
#' @examples
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
#' #Make newdata (3 subjects)
#' newdata <- data.frame(id=rep(1:3, each = 3),x1.1=c(0,0,0,1,0,1,0,1,0),
#'                       x2.2=c(0,1,0,0,0,0,1,0,1), trans = rep(1:3, 3))
#' HvH <- msfit_subjects(cx,newdata,trans=tmat)
#' pt <- probtrans_subjects(HvH,predt=0)
#' 
#' # Default, prediction from state 1
#' summary(pt,id = 1)
#' # Only from states 1 and 3
#' summary(pt,id=2, from=c(1, 3))
#' # Use from=0 for prediction from all states
#' summary(pt,id=3, from=0)
#' # Only to states 1 and 2
#' summary(pt,id=1, to=1:2)
#' # Transition probabilities only at specified time points
#' summary(pt,id=3, times=seq(0, 15, by=3))
#' # Last specified time point is larger than last observed, not printed
#' # Use extend=TRUE as in summary.survfit
#' summary(pt,id=1, times=seq(0, 15, by=3), extend=TRUE)
#' # When the number of time points specified is larger than 12, head and tail is shown
#' x <- summary(pt, id = 2, times=seq(5, 8, by=0.25))
#' x
#' 
#' 
#' @export
summary.probtrans.subjects <- function(object, id, times, from=1, to=0,
                                       extend=FALSE, ...)
{
  if(!inherits(object, "probtrans.subjects")){
    stop("Object to summarize must be of class 'probtrans.subjects'.")
  }
  
  if(length(id) > 1){
    stop("Argument 'id' must be at most of length 1.")
  }
  
  which_id <- which(object[["ids"]] == id)  
  
  if(length(which_id) == 1){
    res <- summary(object[[which_id]], variance = FALSE, extend = extend, times = times,
                   from = from, to = to, ...)
    attr(res, "id") <- id
    return(res)
  } else{
    stop("Subject 'id' not found in 'probtrans.subjects' object. 
         Please provide valid 'id'.")
  }
}
