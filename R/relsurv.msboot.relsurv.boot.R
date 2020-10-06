#' Default theta function used for msboot.relsurv
#' 
#' Helper function used for calling inside msboot.relsurv
#' (used for every bootstrap dataset).
#' This function is used for calculating split hazards
#' and evaluating them at all needed times.
#' @param data An object of class 'msdata' containing a bootstrapped sample
#' @param transmat The transition matrix of class transMat
#' @param all_times All times at which the hazards have to be evaluated
#' @param split.transitions An integer vector containing the numbered transitions that should be split. Use same numbering as in the given transition matrix
#' @param rmap An optional list to be used if the variables in the dataset are not organized (and named) in the same way as in the ratetable object
#' @param time.format Define the time format which is used in the dataset Possible options: c('days', 'years', 'months'). Default is 'days'
#' @param boot_orig_msfit Logical, if true, do the bootstrap for the basic msfit model
#' @param ratetable The population mortality table. A table of event rates, organized as a ratetable object, see for example relsurv::slopop. Default is slopop
#' @param add.times Additional times at which hazards should be evaluated
#' @return A list of calculated values for the given bootstrap sample.
#'
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{msboot.relsurv}}
#' 
#' @export
`msboot.relsurv.boot` <- function(data, transmat, all_times,
                                split.transitions, rmap, time.format, boot_orig_msfit=FALSE, ratetable=relsurv::slopop, add.times){
  if(!missing(rmap))  rmap <- as.call(rmap)
  
  # Round, if needed:
  tolerance <- 1e-15
  # Fit models:
  cox_new <- suppressWarnings(coxph(Surv(Tstart, Tstop, status)~strata(trans), 
                                    data=data, method="breslow"))
  msf_new <- suppressWarnings(mstate::msfit(cox_new,trans=transmat, variance = FALSE))
  msf_relsurv <- msfit.relsurv(msfit.obj = msf_new, data = data,
                               split.transitions = split.transitions
                               ,ratetable = ratetable
                               ,rmap = rmap
                               ,time.format = time.format
                               ,link_trans_ind = TRUE
                               ,substitution = FALSE
                               ,add.times = add.times
  )
  # Hazards:
  output <- msf_relsurv$Haz
  
  output_updated <- output[0,]
  # For every transition in msf_relsurv$Haz, add the times
  # that are missing in this chosen bootstrap sample:
  for(i in msf_relsurv$trans[!is.na(msf_relsurv$trans)]){ 
    # Take transition subsets:
    output_tmp <- subset(output, trans == i)
    all_times_tmp <- subset(all_times, trans == i)
    
    # Find times that have to be added:
    missing_times <- all_times_tmp$Tstop[!(all_times_tmp$Tstop %in% output_tmp$time)]
    missing_times <- unique(sort(missing_times))
    
    # Add missing times, if present at all:
    if(length(missing_times)>0){
      
      # Make new data.frame and add it to the existing one:
      output_tmp_2 <- data.frame(time = missing_times, Haz = NA, trans = i)
      output_tmp <- rbind(output_tmp, output_tmp_2)
      output_tmp <- output_tmp[order(output_tmp$time),]
      
      small_diffs <- c(1, diff(output_tmp$time))
      small_diffs2 <- c(small_diffs[-1], 1)
      
      output_tmp <- output_tmp[!is.na(output_tmp$Haz) | (small_diffs>tolerance & small_diffs2>tolerance),] 
      output_tmp$Haz <- suppressWarnings({mstate:::NAfix(output_tmp$Haz,0)})
    }
    
    # Check if there are any unneeded times:
    unneeded_times <- which(!(output_tmp$time %in% all_times_tmp$Tstop))
    if(length(unneeded_times)>0) output_tmp <- output_tmp[-unneeded_times,]
    
    output_updated <- rbind(output_updated, output_tmp)
    
  }
  if(boot_orig_msfit){
    # ADDITION: BOOTSTRAPPED VAR FROM MSFIT
    output_updated$trans <- paste0(output_updated$trans, "new")
    
    # Find transitions that have been splitted:
    splitted_trans <- msf_new$trans[,!(colnames(msf_new$trans) %in% colnames(msf_relsurv$trans))]
    # Save them in a nice format:
    splitted_trans <- na.omit(as.vector(splitted_trans))
    
    for(i in splitted_trans){
      output_tmp <- subset(msf_new[[1]], trans ==i)
      all_times_tmp <- subset(all_times, trans == i)
      
      missing_times <- all_times_tmp$Tstop[!(all_times_tmp$Tstop %in% output_tmp$time)]
      missing_times <- unique(sort(missing_times))
      
      # Add missing times here:
      if(length(missing_times)>0){
        output_tmp_2 <- data.frame(time = missing_times, Haz = NA, trans = i)
        
        output_tmp <- rbind(output_tmp, output_tmp_2)
        output_tmp <- output_tmp[order(output_tmp$time),]
        
        small_diffs <- c(1, diff(output_tmp$time))
        small_diffs2 <- c(small_diffs[-1], 1)
        
        output_tmp <- output_tmp[!is.na(output_tmp$Haz) | (small_diffs>tolerance & small_diffs2>tolerance),]
        # output_tmp <- output_tmp[small_diffs>tolerance,]
        
        output_tmp$Haz <- mstate:::NAfix(output_tmp$Haz,0)
      }
      output_tmp$trans <- paste0(output_tmp$trans, "old")
      output_updated <- rbind(output_updated, output_tmp)
    }
  }
  return(list(output_updated, msf_relsurv$link_trans))
}