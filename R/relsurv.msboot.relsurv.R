#' Bootstrap function for upgraded multi-state models using relsurv
#' 
#' A helper nonparametric bootstrapping function for variances
#' in extended multi-state models using relative survival.
#' This implementation is written based on function mstate:::msboot.
#' @param theta	A function of data and perhaps other arguments, returning the value of the statistic to be bootstrapped
#' @param data An object of class 'msdata', such as output from msprep
#' @param B	The number of bootstrap replications; the default is taken to be quite small (5) since bootstrapping can be time-consuming
#' @param id Character string indicating which column identifies the subjects to be resampled
#' @param verbose The level of output; default 0 = no output, 1 = print the replication
#' @param transmat The transition matrix of class transMat
#' @param all_times All times at which the hazards have to be evaluated
#' @param split.transitions An integer vector containing the numbered transitions that should be split. Use same numbering as in the given transition matrix
#' @param rmap An optional list to be used if the variables in the dataset are not organized (and named) in the same way as in the ratetable object
#' @param time.format Define the time format which is used in the dataset Possible options: c('days', 'years', 'months'). Default is 'days'
#' @param boot_orig_msfit Logical, if true, do the bootstrap for the basic msfit model
#' @param ratetable The population mortality table. A table of event rates, organized as a ratetable object, see for example relsurv::slopop. Default is slopop
#' @param add.times Additional times at which hazards should be evaluated
#' @param ... Any further arguments to the function theta
#' @return A list of size B containing the results for every bootstrap iteration.
#' 
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}, Marta Fiocco, Hein Putter \email{H.Putter@@lumc.nl}
#' @seealso \code{\link{msboot}}
#' 
#' @export
`msboot.relsurv` <- function(theta, data, B = 5, id = "id", verbose = 0, 
                           transmat, all_times, split.transitions,
                           rmap, time.format, boot_orig_msfit, ratetable=relsurv::slopop, add.times, ...){
  
  if (!inherits(data, "msdata")) 
    stop("'data' must be a 'msdata' object")
  trans <- attr(data, "trans")
  ids <- unique(data[, id])
  n <- length(ids)
  th <- theta(data, transmat=transmat, all_times=all_times, split.transitions=split.transitions,
              rmap=rmap, time.format=time.format, boot_orig_msfit=boot_orig_msfit, ratetable=ratetable, add.times=add.times, ...)
  nr <- nrow(th)
  # Prepare res object:
  res <- vector("list", B)
  
  find_trans <- as.numeric(stats::na.omit(as.vector(transmat)))
  take_values <- rep(TRUE, B)
  
  add.times.ind <- FALSE
  # Take care of add.times:
  if(!missing(add.times)){
    add.times.ind <- TRUE
  }
  
  # start_it <- 1
  for (b in 1:B) {
    if (verbose > 0) {
      cat("\nBootstrap replication", b, "\n")
      flush.console()
    }
    bootdata <- NULL
    bids <- sample(ids, replace = TRUE)
    bidxs <- unlist(sapply(bids, function(x) which(x == 
                                                     data[, id])))
    bootdata <- data[bidxs, ]
    if (verbose > 0) {
      print(date())
      print(events(bootdata))
      cat("applying theta ...")
    }
    
    if(length(unique(bootdata[,id]))==1){
      take_values[b] <- FALSE
      next
    } 
    if(add.times.ind){
      add.times1 = add.times
      add.times1 = add.times1[add.times1 <= max(bootdata$Tstop, na.rm=TRUE)]
      thstar <- theta(data=bootdata, transmat=transmat, all_times=all_times, 
                      split.transitions=split.transitions,
                      rmap=rmap, time.format=time.format, boot_orig_msfit=boot_orig_msfit, ratetable=ratetable, add.times=add.times1, ...)
    } else{
      thstar <- theta(data=bootdata, transmat=transmat, all_times=all_times, 
                      split.transitions=split.transitions,
                      rmap=rmap, time.format=time.format, boot_orig_msfit=boot_orig_msfit, ratetable=ratetable, add.times=add.times, ...)
    }
    
    # Check if some transitions have to be removed 
    # (if nobody can make that transition, i.e. if there's no one at the at-risk set)
    observed_trans <- sort(unique(bootdata$trans))
    if(!identical(observed_trans,find_trans)){
      which_trans <- find_trans[!(find_trans %in% observed_trans)]
      new_trans <- unlist(thstar[[2]][which_trans])
      thstar[[1]] <- subset(thstar[[1]], !(trans %in% new_trans))
    }
    
    thstar <- thstar[[1]]
    thstar$b <- b
    
    # Save result: 
    res[[b]] <- thstar
  }
  if (verbose) 
    cat("\n")
  
  res <- res[take_values]
  return(res)
}
