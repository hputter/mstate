#' Compute subject-specific transition probabilities
#' 
#' This function computes subject-specific transition probabilities
#' in multi-state models.
#' 
#' For details refer to de Wreede, Fiocco & Putter (2010).
#' 
#' @param object \link[=msfit_subjects]{msfit.subjects} object containing estimated cumulative hazards
#' for each of the transitions in the multi-state model
#' @param predt A positive number indicating the prediction time. This is
#' either the time at which the prediction is made (if \code{direction}=
#' \code{"forward"}) or the time for which the prediction is to be made (if
#' \code{direction}=\code{"fixedhorizon"})
#' @param direction One of \code{"forward"} (default) or \code{"fixedhorizon"},
#' indicating whether prediction is forward or for a fixed horizon
#'
#'
#' @return An object of class \code{"probtrans.subjects"}.
#' This is a list of length n (number of subjects in \code{object$ids}), with each list element
#' an object of class \code{\link{probtrans}} for the associated 
#' subject. List elements can be accessed using \code{[[x]]}, with \code{x} 
#' ranging from 1 to n. Additionally, each list element 
#' has an element \code{$id}, representing the subject id and the output object 
#' also has an element \code{$ids} representing the subject ids in order. 
#' Plot and summary methods have been defined for
#' \code{"probtrans(.subjects)"} objects.
#' 
#' 
#' @author Daniel Gomon
#' @references Andersen PK, Borgan O, Gill RD, Keiding N (1993).
#' \emph{Statistical Models Based on Counting Processes}. Springer, New York.
#' 
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: Competing
#' risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
#' 2389--2430.
#' 
#' Therneau TM, Grambsch PM (2000). \emph{Modeling Survival Data: Extending the
#' Cox Model}. Springer, New York.
#' 
#' de Wreede LC, Fiocco M, and Putter H (2010). The mstate package for
#' estimation and prediction in non- and semi-parametric multi-state and
#' competing risks models. \emph{Computer Methods and Programs in Biomedicine}
#' \bold{99}, 261--274.
#' 
#' de Wreede LC, Fiocco M, and Putter H (2011). mstate: An R Package for the
#' Analysis of Competing Risks and Multi-State Models. \emph{Journal of
#' Statistical Software}, Volume 38, Issue 7.
#' @keywords survival
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
#'         data=tglong,method="breslow")
#' summary(cx)
#' #Make newdata (3 subjects)
#' newdata <- data.frame(id=rep(1:3, each = 3),x1.1=c(0,0,0,1,0,1,0,1,0),
#'                       x2.2=c(0,1,0,0,0,0,1,0,1), trans = rep(1:3, 3))
#' HvH <- msfit_subjects(cx,newdata,trans=tmat)
#' # probtrans_subjects
#' pt <- probtrans_subjects(HvH,predt=0)
#' # A 'probtrans' object for person 1
#' pt[[1]]
#' 
#' @export 
`probtrans_subjects` <- function(object,predt,direction=c("forward","fixedhorizon")){
  
  #Check if input is a msfit.subjects object
  if(!inherits(object, "msfit.subjects")){
    stop("object must be an object of class 'msfit.subjects'.")
  }
  
  #Match direction
  direction <- match.arg(direction, choices = c("forward", "fixedhorizon"))
  
  #Prediction time must be positive!
  if(!(predt >= 0) | !(is.numeric(predt))){
    stop("Prediction time must be a positive numeric value.")
  }
  
  #--------Step 4: Calculate transition probabilities per subject------#
  
  #This step is relatively easy. We simply use probtrans_D for each subject separately.
  
  n_subjects <- length(object$ids)
  
  
  #We now determine n_subjects and subject_ids in the beginning of this function.
  res <- vector(mode = "list", length = n_subjects)
  for(i in 1:n_subjects){
    res[[i]] <- probtrans_D(list(intensity_matrices = object$intensities[, , , i],
                                 unique_times = object$times), 
                            predt = predt, direction = direction, as.df = TRUE)
    class(res[[i]]) <- "probtrans"
    res[[i]]$trans <- object$trans
    res[[i]]$method <- "aalen"
    res[[i]]$predt <- predt
    res[[i]]$direction <- direction
    res[[i]]$id <- object$ids[i]
  }
  res$ids <- object$ids
  class(res) <- "probtrans.subjects"
  #-----------OUTPUT------------#
  return(res)
}