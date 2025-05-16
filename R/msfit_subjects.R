#' Compute subject-specific transition hazards for multiple subjects at once
#' 
#' This function computes subject-specific or overall cumulative transition
#' hazards for each of the possible transitions in the multi-state model.
#' Contrary to \code{\link[msfit]{msfit}}, this function allows to 
#' calculate hazards for multiple subjects at once, but does not allow
#' to calculate (co)variances (yet).
#' 
#' @details
#' The data frame needs to have one row for each transition in the multi-state
#' model, per subject: (n_subjects x n_transitions) rows in total.
#' Contrary to \code{\link[msfit]{msfit}}, it is not necessary to 
#' manually specify the \code{strata} variable,
#' as long as the strata can be determined from the data using the formula 
#' used in \code{object}. For details refer to de Wreede, Fiocco & Putter (2010). 
#' So far,
#' the results have been checked only for the \code{"breslow"} method of
#' dealing with ties in \code{\link[survival:coxph]{coxph}}, so this is
#' recommended.
#' The goal of this function is to determine the subject specific intensity matrices
#' \eqn{dA_i(t)}{dA_i(t)}, which can be further used to determine transition probabilities
#' through the product integral \deqn{P_i(s,t) = \prod_{(s,t]}(I + dA_i(u))}{P_i(s,t) = \prod_{(s,t])}(I + dA_i(u))}
#' 
#' 
#' 
#' 
#' @param object A \code{\link[survival:coxph]{coxph}} object describing the
#' fit of the multi-state model, must contain a 'strata()' term.
#' @param newdata A \code{data.frame} containing K rows per
#' subject in the data, with K the number of transitions possible in the model.
#' The rows pertaining to a single subject must be ordered according to the 
#' transition numbers in \code{trans}.
#' The following named columns must be present:
#' \describe{
#'   \item{\code{id}:}{Unique identifier of the subject, must be numeric or character;}
#'   \item{\code{"variables"}:}{The subject-specific covariates (appearing in the \code{coxph} formula),
#'   ordered according to transition number in \code{trans}.}
#' } Note that newdata must contain a column containing the variable which was 
#' used to determine the stratum of a transition in \code{object}. 
#' @param trans Transition matrix describing the states and transitions in the
#' multi-state model. See \code{trans} in \code{\link{msprep}} for more
#' detailed information
#' 
#' 
#' @return An object of class \code{"msfit.subjects"}, which is a list containing
#' \item{intensities}{A 4-dimensional \code{array} containing the estimated cumulative 
#' intensity increments at \code{times}. The first and second dimensions pertain 
#' to the to/from states of the transitions, the third to the \code{times} and 
#' the fourth to the unique subject \code{ids}.} 
#' \item{times}{The unique times at which the intensity increments were determined.}
#' \item{ids}{The unique id's used to represent the individual subjects, taken 
#' from the 'id' column in \code{newdata}.} 
#' \item{trans}{The transition matrix used.}
#' Looking for an 'msfit' object for a single subject? See 
#' \code{\link{msfit_subjects_to_msfit}}.
#' 
#' 
#' @author Daniel Gomon
#' 
#' @seealso \code{\link{plot.msfit}}, \code{\link{msfit_subjects_to_msfit}}
#' 
#' @references Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
#' Competing risks and multi-state models. \emph{Statistics in Medicine}
#' \bold{26}, 2389--2430.
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
#' # expanded covariates
#' tglong <- expand.covs(tglong,c("x1","x2"))
#' # Cox model with different covariate
#' cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
#' 	data=tglong,method="breslow")
#' #Fit on multiple subjects at once, by providing 'id' column.
#' #cx was fit using strata(trans), so we must have a 'trans' column.
#' newdata <- data.frame(id=rep(1:3, each = 3),x1.1=c(0,0,0,1,0,1,0,1,0),
#'                       x2.2=c(0,1,0,0,0,0,1,0,1), trans = rep(1:3, 3))
#' msf_subj <- msfit_subjects(cx,newdata,trans=tmat)
#' #Extract an 'msfit' object for subject 1
#' msf_subj1 <- msfit_subjects_to_msfit(msf_subj, 1)
#' #Can now use standard 'msfit' plotting capabilities
#' plot(msf_subj1)
#' 
#' @export 
msfit_subjects <- function(object, newdata, trans){
  
  
  #---------------DATA CHECKS----------------#
  
  #We first copy some checks from msfit, which are actually from survfit
  #Here we make sure that the cox fit is appropriate, and extract the 
  #correct terms for each transition.
  
  #Some checks on coxmod
  if(!is.null((object$call)$weights) || !is.null(object$weights))
    stop("msfit cannot (yet) compute the result for a weighted model")
  Terms <- terms(object)
  strat <- attr(Terms, "specials")$strata
  if (is.null(strat)) stop("object has to have strata() term")
  cluster <- attr(Terms, "specials")$cluster
  if (length(cluster)) stop("cluster terms are not supported")
  if (!is.null(attr(object$terms, "specials")$tt))
    stop("msfit cannot yet process coxph models with a tt term")
  
    
  if(!inherits(trans, "matrix")){
    stop("trans must be a transition matrix (array).")
  }
  
  if(!is.data.frame(newdata)){
    stop("newdata must be a data frame containing 'id' and variable columns, one row per transition.")
  }
  
  tmat2 <- to.trans2(trans)[, c(2,3,1)]
  names(tmat2)[3] <- "trans"
  n_transitions <- nrow(tmat2)
  
  
  #If no "id" column is supplied, assume all entries are from the same subject.
  #If so, there need to be exactly n_transitions rows in the data
  if(!"id" %in% colnames(newdata)){
    if(nrow(newdata) == n_transitions){
      newdata <- cbind(newdata, data.frame(id = rep(1, n_transitions)))
      warning("No 'id' column specified, assuming all rows relate to a single subject.")
    } else{
      stop("Please identify subject rows in 'newdata' using a column named 'id'.")
    }
  }
  
  #Check whether id is character or numeric or integer
  if(!inherits(newdata[, "id"], "numeric") && !inherits(newdata[, "id"], "character") 
     && !inherits(newdata[, "id"], "integer")){
    stop("The 'id' column in 'newdata' must be of class numeric (or integer) or character.")
  }
  
  #Order according to 'id'
  newdata <- newdata[order(newdata[, "id"]), ]
  
  #Keep track of id's and number of subjects
  subject_ids <- unique(newdata[, "id"])
  n_subjects <- length(subject_ids)
  

  
  
  
  
  #----------Step 1: Baseline intensities-------------#
  #Recover baseline intensities from cox model
  #Output from get_intensity matrices
  #List with $intensity_matrices and $unique_times
  #Automatically determines strata used in the coxmodel from the formula.
  #Checks to see if the named column in the cox "object" is present in newdata as well
  baseline_intensities <- baseline_intensities_from_coxmod(object, trans)
  
  #---------Step 2: Subject specific risks-------------#
  
  #####Step 2.2: Extract the subject specific risk from the transformed data#####
  
  subject_specific_risks <- trans_specific_risks(object = object, newdata = newdata,
                                                 trans = trans, n_subjects = n_subjects,
                                                 subject_ids = subject_ids)

  #---------Step 3: Obtain subject specific predictions----------#
  
  #To obtain subject specific predictions we simply multiply the risk matrix
  #with the appropriate indices of the baseline intensity matrices
  #print("Step 3")
  subject_specific_intensity_matrices <- subject_specific_intensity_matrices(subject_specific_risks = subject_specific_risks,
                                                                             baseline_intensities = baseline_intensities,
                                                                             trans = trans)
  
  out <- list(intensities = subject_specific_intensity_matrices$subject_intensity_matrices,
              times = subject_specific_intensity_matrices$unique_times,
              ids = subject_ids,
              trans = trans)
  class(out) <- "msfit.subjects"
  return(out)
}


#' Obtain an 'msfit' object for a single subject from an 'msfit.subjects' 
#' object.
#' 
#' @param object A \code{'msfit.subjects'} object.
#' @param id A single subject_id, must be contained in \code{object$ids}.
#' 
#' @returns An \code{\link[mstate:msfit]{msfit}} object for a single subject.
#' 
#' @seealso \code{\link{msfit}}
#' 
#' 
#' @author Daniel Gomon
#' 
#' 
#' 
#' @export
#' 

msfit_subjects_to_msfit <- function(object, id){
  
  #Check class of object
  if(!inherits(object, "msfit.subjects")){
    stop("object must be an 'msfit.subjects' object.")
  }
  #Check if valid id
  if(!id %in% object$ids){
    stop("id must be contained in 'object$ids'.")
  }
  
  #Find position of id within intensity matrices:
  subj_pos <- which(object$ids == id)
  
  #Get tmat2
  tmat2 <- to.trans2(object$trans)
  n_transitions <- nrow(tmat2)
  n_states <- nrow(object$trans)
  
  #Translate subject-specific intensity matrix to msfit model.
  times <- object$times
  
  #Extract all cumulative hazards (simple cumsum() application on 3rd dimension 
  #of subject specific intensity matrix)
  Haz <- apply(apply(object$intensities[, , , subj_pos], 3, 
                     function(x) x), 1, cumsum)
  #Retain only relevant column positions
  relevant_columns <- (tmat2$to-1)*n_states + tmat2$from
  Haz <- c(Haz[, relevant_columns])
  
  cHaz <- data.frame(time = rep(times, n_transitions),
                     Haz = Haz,
                     trans = rep(tmat2$transno, each = length(object$times)))
  
  #Initialize msfit object
  out <- list(Haz = cHaz,
              trans = object$trans)
  class(out) <- "msfit"
  out
}




# Internal Functions of msfit_subjects ------------------------------------



#' Recover baseline intensities (in the form of intensity_matrices) from a 
#' coxph() fit on multi-state data.
#' 
#' @inherit msfit_subjects params
#' @param tmat A transition matrix as created by \code{\link[mstate:transMat]{transMat}}.
#' 
#' @importFrom mstate msfit to.trans2
#' @import survival
#' @importFrom stats formula model.frame as.formula predict
#' 
#' @author Daniel Gomon
#' 
#' @keywords internal


baseline_intensities_from_coxmod <- function(object, tmat){
  #Extract baseline intensities from a cox model
  #Before running this function, check whether the cox model has been fit 
  #for a MSM (contains strata & no interaction/tt terms)
  
  #object must be a coxph model with strata
  #tmat must be a transition matrix
  #Create a baseline subject, which has 0 variables everywhere!
  #We need strata here, as we want the baseline hazard for each strata separately
  ttmat <- to.trans2(tmat)[, c(2, 3, 1)]
  #Setting the name to "trans" is in line with msprep... This kind of requires people to use msprep.
  names(ttmat)[3] <- "trans"
  baseline_subject <- matrix(0, ncol = length(object$coefficients) + 3, 
                             nrow = nrow(ttmat), 
                             dimnames = list(NULL, c(colnames(ttmat), 
                                                     names(object$coefficients))))
  baseline_subject <- as.data.frame(baseline_subject)
  baseline_subject[1:ncol(ttmat), 1:nrow(ttmat)] <- ttmat
  #Add strata to the subject observations
  #Check how the strata were created
  Terms <- terms(object)
  stangle <- untangle.specials(Terms, 'strata')
  #stangle$vars is now the formula for strata creation
  strata_column <- model.frame(as.formula(paste0(stangle$vars, " ~ 1")), 
                               data = baseline_subject)
  baseline_subject <- cbind(baseline_subject, strata_column)
  
  #Use msfit once to obtain baseline hazards
  msfit_baseline <- mstate::msfit(object, newdata = baseline_subject, trans = tmat)
  
  #Now extract the intensities into matrix form
  baseline_intensities <- get_intensity_matrices(msfit_baseline)
  
  
  
  out <- baseline_intensities
  return(out)
}



#' Calculate subject specific risks for subjects in newdata
#' 
#' @description Return a 2 dimensional array, with subjects in the first dimension and
#' transition (numbers) in the second. Can expand later to include time-dependent
#' covariates by introducing extra dimension. Entries of the matrix are exp(lp)_i^m,
#' with i denoting the subject and m the transition.
#' 
#' @inherit msfit_subjects params
#' 
#' @import survival
#'
#' @author Daniel Gomon
#' 
#' @keywords internal
#' 

trans_specific_risks <- function(object, newdata, trans, n_subjects, subject_ids){
  #Given a coxph object (fitted on multi-state data)
  #We want to extract the subject-specific risk (exp(lp)_1, exp(lp)_2, ..., exp(lp)_k)
  #for each transition and each unique subject in the data.
  
  #object = coxph() fit
  #newdata = data.frame/matrix with covariates IN EXTENDED FORM. Must contain id 
  #trans = transition matrix
  #strata_column = a column vector indicating which transitions belong to which strata,
  #must match with the to.trans2(trans) ordering.
  
  #Required output:
  #2D array with:
  #rows: subjects i = 1, ..., n
  #columns: transitions m = 1, ..., M
  #3rd dimension (future): unique times t = t_1, ..., t_K
  #Each entry of the array should be exp(lp)_{i, m, t}
  #so the risk for person i in transition m (at time t (extension, consider later))
  
  
  
  tmat2 <- to.trans2(trans)
  n_transitions <- nrow(tmat2)
  
  
  #Initialize output matrix
  subject_specific_risks_mat <- matrix(NA, nrow = n_subjects, ncol = n_transitions)
  rownames(subject_specific_risks_mat) <- subject_ids
  
  #Determine risk for each subject separately
  for(i in seq_along(subject_ids)){
    tempdat <- newdata[newdata[, "id"] == subject_ids[i],]
    #For each transition, we now obtain the subject risk in a vector of length n_transitions
    risk_subject <- predict(object = object, newdata = tempdat, type = "risk",
                            reference = "zero")
    subject_specific_risks_mat[i,] <- risk_subject
  }
  
  return(subject_specific_risks_mat)
}


#' Calculate the subject specific intensity matrices
#' 
#' @description
#' For each subject, calculate a 3D array containing the states (from, to) in the
#' first two dimension and the times in the third.
#' This is then stored in a 4D array, with the 4th dimension indicating each unique subject.
#' 
#' 
#' @author Daniel Gomon
#' 
#' 
#' @keywords internal

subject_specific_intensity_matrices <- function(subject_specific_risks, baseline_intensities, trans){
  #Given subject specific risks, we want to obtain the intensity matrices for each subject
  
  
  tmat2 <- to.trans2(trans)
  n_subjects <- nrow(subject_specific_risks)
  n_transitions <- nrow(tmat2)
  n_times <- length(baseline_intensities$unique_times)
  
  #pre-specify the output
  #First 3 dimensions: from, to, time
  #Fourth dimension: subjects
  intensity_matrices_zero_diag <- baseline_intensities$intensity_matrices
  for(k in 1:n_times){
    diag(intensity_matrices_zero_diag[, , k]) <- 0
  }
  subject_intensity_matrices <- array(intensity_matrices_zero_diag, 
                                      dim = c(dim(intensity_matrices_zero_diag), n_subjects))
  dimnames(subject_intensity_matrices)[[4]] <- rownames(subject_specific_risks)
  
  for(i in 1:n_subjects){
    for(m in 1:n_transitions){
      from <- tmat2$from[m]
      to <- tmat2$to[m]
      subject_intensity_matrices[from, to, , i] <- subject_intensity_matrices[from, to, , i] * subject_specific_risks[i, m]
    }
    for(k in 1:n_times){
      diag(subject_intensity_matrices[, , k, i]) <- 1 - rowSums(subject_intensity_matrices[, , k, i])
    }
  }
  
  out <- list(subject_intensity_matrices = subject_intensity_matrices,
              unique_times = baseline_intensities$unique_times)
  
  return(out)
  
}

