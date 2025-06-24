#' Plot method for an msfit.subjects object
#' 
#' Plot method for an object of class \code{"msfit.subjects"}. It plots the estimated
#' cumulative transition intensities in the multi-state model for a single subject.
#' Wrapper for \code{\link{plot.msfit}}.
#' 
#' @param x Object of class \code{"msfit.subjects"}, containing estimated transition intensity matrices
#' for all transitions in a multi-state model for (multiple) subjects
#' @param id Identifier of subject to make plot for. 
#' @param \dots Further arguments to \code{\link{plot.msfit}}
#' 
#' @method plot msfit.subjects
#' 
#' @return No return value
#' 
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#' @author Daniel Gomon
#'
#' @seealso \code{\link{msfit_subjects}}
#' @keywords hplot
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # data in wide format, for transition 1 this is dataset E1 of
#' # Therneau & Grambsch (2000)
#' tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'         dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'         x1=c(1,1,1,0,0,0),x2=c(6:1))
#' # data in long format using msprep
#' tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
#' 		data=tg,keep=c("x1","x2"),trans=tmat)
#' # events
#' events(tglong)
#' table(tglong$status,tglong$to,tglong$from)
#' # expanded covariates
#' tglong <- expand.covs(tglong,c("x1","x2"))
#' # Cox model with different covariate
#' cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
#' 	data=tglong,method="breslow")
#' summary(cx)
#' #Fit on multiple subjects at once, by providing 'id' column.
#' #cx was fit using strata(trans), so we must have a 'trans' column.
#' newdata <- data.frame(id=rep(1:3, each = 3),x1.1=c(0,0,0,1,0,1,0,1,0),
#'                       x2.2=c(0,1,0,0,0,0,1,0,1), trans = rep(1:3, 3))
#' msf_subj <- msfit_subjects(cx,newdata,trans=tmat)
#' # standard plot for subject 1
#' plot(msf_subj,id = 1)
#' # standard plot for subject 2
#' plot(msf_subj,id = 2)
#' # specifying line width, color, and legend
#' plot(msf_subj,id = 1,lwd=2,col=c("darkgreen","darkblue","darkred"),legend=c("1->2","1->3","2->3"))
#' # separate plots for each transition
#' par(mfrow=c(2,2))
#' plot(msf_subj, id = 1, type="separate",lwd=2)
#' par(mfrow=c(1,1))
#' 
#' @export 
plot.msfit.subjects <- function(x, 
                       id,
                       ...) {
  
  # Prelim
  if (!inherits(x, "msfit.subjects")) stop("'x' must be a 'msfit.subjects' object")
  
  #Transform to msfit object
  x <- msfit_subjects_to_msfit(x, id = id)
  
  #Use existing plotting functionality
  plot(x, ...)
}
