#' @title
#' Plot an object of class \code{"probtrans.subjects"}
#' 
#' @description
#' Plots the transition probabilities for a specific subject. Wrapper for 
#' \code{\link[mstate:plot.probtrans]{plot.probtrans}}
#' 
#' @param x An object of class \code{"probtrans.subjects"}
#' @param id Subject identifier
#' @param ... Further arguments to \code{\link[mstate:plot.probtrans]{plot.probtrans}}
#' 
#' @details
#' Note that 
#'
#' @author Hein Putter and Daniel Gomon
#' 
#' @import mstate
#' 
#' @method plot probtrans.subjects
#' @export
#' 
#' @examples
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # data in wide format, for transition 1 this is dataset E1 of
#' # Therneau and Grambsch (2000)
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
#' # new data, to check whether results are the same for transition 1 as
#' # those in appendix E.1 of Therneau and Grambsch (2000)
#' newdata <- data.frame(id=rep(1:3, each = 3),x1.1=c(0,0,0,1,0,1,0,1,0),
#'                       x2.2=c(0,1,0,0,0,0,1,0,1), trans = rep(1:3, 3))
#' msf <- msfit_subjects(cx,newdata,trans=tmat)
#' # probtrans
#' pt <- probtrans_subjects(msf,predt=0)
#' # default plot
#' plot(pt,id=2,ord=c(2,3,1),lwd=2,cex=0.75)
#' # filled plot
#' plot(pt,id=3,type="filled",ord=c(2,3,1),lwd=2,cex=0.75)
#' # single plot
#' plot(pt,id=1,type="single",lwd=2,col=rep(1,3),lty=1:3,legend.pos=c(8,1))
#' # separate plots
#' par(mfrow=c(2,2))
#' plot(pt,id=1,type="sep",lwd=2)
#' par(mfrow=c(1,1))
#' 
#' # ggplot version - see vignette for details
#' library(ggplot2)
#' plot(pt,id=1,ord=c(2,3,1), use.ggplot = TRUE)
#' 


plot.probtrans.subjects <- function(x, id, ...){
  if(!inherits(x, "probtrans.subjects")){
    stop("Object to plot must be of class 'probtrans.subjects'.")
  }
  
  if(length(id) > 1){
    stop("Argument 'id' must be at most of length 1.")
  }
  
  which_id <- which(x[["ids"]] == id)  
  
  if(length(which_id) == 1){
    plot(x[[which_id]], ...)
  } else{
    stop("Subject 'id' not found in 'probtrans.subjects' object. 
         Please provide valid 'id'.")
  }
}