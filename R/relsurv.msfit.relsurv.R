#' Extend a multi-state model using relative survival 
#' 
#' A function that takes a fitted msfit object and upgrades
#' it using relative survival, where chosen transitions are
#' split in population and excess transitions. This upgraded 
#' msfit object contains the split hazards based on the transition
#' matrix (transMat). The (co)variance matrix is also upgraded, if provided.
#' @param msfit.obj The msfit object which has to be upgraded
#' @param data The data used for fitting the msfit model
#' @param split.transitions An integer vector containing the numbered transitions that should be split. Use same numbering as in the given transition matrix
#' @param ratetable The population mortality table. A table of event rates, organized as a ratetable object, see for example relsurv::slopop. Default is slopop
#' @param rmap An optional list to be used if the variables in the data are not organized (and named) in the same way as in the ratetable object
#' @param time.format Define the time format which is used in the data. Possible options: c('days', 'years', 'months'). Default is 'days'
#' @param var.pop.haz If 'fixed' (default), the Greenwood estimator for the variances is used, where it is assumed that the variance of the population hazards is zero. If 'bootstrap', one gets boostrap estimates for all all transitions. Option 'both' gives both variance estimates
#' @param B Number of bootstrap replications. Relevant only if var.pop.haz == 'bootstrap' or 'both'. Default is B=10.
#' @param seed Set seed
#' @param add.times Additional times at which hazards should be evaluated
#' @param substitution Whether function substitute should be used for rmap argument. Default is TRUE
#' @param link_trans_ind Choose whether the linkage between the old and new transition matrix should be saved. Default is FALSE.
#' @return Returns a msfit object that contains estimates for the extended model
#' with split (population and excess) transitions.
#' 
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{msfit}}
#' @references Manevski D, Putter H, Pohar Perme M, Bonneville EF, Schetelig J, de Wreede LC (2021).
#' Integrating relative survival in multi-state models -- a non-parametric approach.
#' https://arxiv.org/abs/2106.12399
#' 
#' @examples 
#' 
#' \dontrun{
#' library(mstate)
#' # Load dataset:
#' data("ebmt1")
#' # Transition matrix:
#' tmat <- transMat(list(c(2,3),c(4), c(), c()), 
#'                  names = c("Alive relapse-free", "Relapse","NRM", "DaR"))
#' # Data in long format using msprep
#' df <- msprep(time=c(NA,"rel","srv","srv"), status=c(NA,"relstat","srvstat","srvstat"),
#'              data=ebmt1, trans=tmat)
#' # Generate demographic covariates (which are usually present in datasets) 
#' # and based on them estimate the population hazard.
#' set.seed(510)
#' df$age <- runif(nrow(df), 45, 65)
#' df$sex <- sample(c("male", "female"), size = nrow(df), replace = TRUE)
#' df$dateHCT <- sample(seq(as.Date('1990/01/01'), 
#'     as.Date('2000/01/01'), by="day"), nrow(df), replace = TRUE) # generate years
#' # Cox object:
#' cx <- coxph(Surv(Tstart,Tstop,status)~strata(trans),
#'             data=df,method="breslow")
#' # Basic multi-state model:
#' mod <- msfit(cx,trans=tmat)
#' # Extended multi-state model, where the two transition
#' # reaching death are split in excess and population parts.
#' # We assume patients live like in the Slovene population,
#' # thus we use Slovene mortality tables in this example.
#' # Variances estimated using 25 bootstrap replications.
#' mod.relsurv <- msfit.relsurv(msfit.obj = mod, data=df, split.transitions = c(2,3),
#'                             ratetable = relsurv::slopop, 
#'                             rmap = list(age=age*365.241, year=dateHCT),
#'                             time.format = "days",
#'                             var.pop.haz = "bootstrap",
#'                             B = 25)
#' # Estimate transition probabilities:
#' pt <- probtrans(mod.relsurv, predt=0, method='greenwood')
#' # Estimated cumulative hazards with the corresponding 
#' # bootstrap standard errors at 300, 600, 900 days:
#' summary(object = mod.relsurv, times = c(300, 600, 900), conf.type = 'log')
#' # Estimated transition probabilities together with the corresponding 
#' # bootstrap standard errors and log.boot confidence intervals 
#' # at 300, 600, 900 days:
#' summary(object = pt, times = c(300, 600, 900), conf.type = 'log')
#' # Plot the measures:
#' plot(mod.relsurv, use.ggplot = TRUE)
#' plot(pt, use.ggplot = TRUE)
#' }
#' @export
`msfit.relsurv` <- function(msfit.obj, 
                          data,
                          split.transitions, 
                          ratetable = relsurv::slopop,
                          rmap, 
                          time.format = "days",
                          var.pop.haz = c("fixed", "bootstrap", "both"),
                          B = 10,
                          seed = NULL,
                          add.times,
                          substitution=TRUE,
                          link_trans_ind = FALSE
){
  # NOTE: var.pop.haz="fixed" gives covariances. Bootstrap does not.
  
  if (!requireNamespace("relsurv", quietly = TRUE)) {
    stop("Package \"relsurv\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  ####################
  # 1. Prepare objects
  ####################
  trans1 <- NULL
  trans2 <- NULL
  b <- NULL
  
  if(!missing(rmap)){
    if(substitution){
      rmap <- substitute(rmap)
    }
  }
  
  var.pop.haz <- match.arg(var.pop.haz)
  bootstrap_ind <- FALSE # helper indicator
  
  # Add additional times, if needed:
  if(!missing(add.times)){
    if(max(add.times) > max(msfit.obj[[1]]$time)){
      warning("Maximum value in add.times is bigger than maximum time in dataset. Do you really want to extrapolate the hazards? \n")
    }
    
    add.times <- add.times[!(add.times %in% msfit.obj[[1]]$time)]
    
    msfit.new <- msfit.obj[[1]][0,]
    # Make new data.frame and add it to the existing one:
    for(tr in unique(msfit.obj[[1]]$trans)){
      msfit_tmp <- subset(msfit.obj[[1]], trans==tr)
      
      msfit_tmp_2 <- data.frame(time = add.times, Haz = NA, trans = tr)
      msfit_tmp <- rbind(msfit_tmp, msfit_tmp_2)
      msfit_tmp <- msfit_tmp[order(msfit_tmp$time),]
      
      msfit_tmp$Haz <- NAfix(msfit_tmp$Haz,0)
      
      msfit.new <- rbind(msfit.new, msfit_tmp)
    }
    msfit.obj[[1]] <- msfit.new
    
    # Take care of varHaz, if needed:
    if(length(msfit.obj) == 3){
      msfit.new <- msfit.obj[[2]][0,]
      trans.grid <- unique(msfit.obj[[2]][,c("trans1", "trans2")])
      # Make new data.frame and add it to the existing one:
      for(tr in 1:nrow(trans.grid)){
        msfit_tmp <- subset(msfit.obj[[2]], (trans1==trans.grid[tr, "trans1"] & trans2==trans.grid[tr, "trans2"]))
        
        msfit_tmp_2 <- data.frame(time = add.times, varHaz = NA, trans1 = trans.grid[tr, "trans1"], trans2 = trans.grid[tr, "trans2"])
        msfit_tmp <- rbind(msfit_tmp, msfit_tmp_2)
        msfit_tmp <- msfit_tmp[order(msfit_tmp$time),]
        
        msfit_tmp$varHaz <- NAfix(msfit_tmp$varHaz,0)
        
        msfit.new <- rbind(msfit.new, msfit_tmp)
      }
      msfit.obj[[2]] <- msfit.new
    }
  }
  
  # Exctract hazards:
  Haz <- msfit.obj[[1]]
  # Exctract covariances, if present:
  if(length(msfit.obj) == 3){
    varHaz <- msfit.obj[[2]]
  }
  
  # Define time-related objects:
  Year <- 365.241
  Month <- Year/12
  
  data_orig <- data
  time.format.orig <- time.format
  
  if(time.format == "days"){
    if(max(msfit.obj[[1]]$time) < 30) warning("Your max time in the data is less than 30 days. If time is not stored in days, please use argument time.format. \n")
  }
  else if(time.format == "years"){
    if(max(msfit.obj[[1]]$time) > 100) warning("Your max time in the data is more than 100 years. If time is not stored in years, please use argument time.format. \n")
    
    data$Tstart <- data$Tstart*Year
    data$Tstop <- data$Tstop*Year
    data$time <- data$time*Year
    
    Haz$time <- Haz$time*Year
    if(exists("varHaz")) varHaz$time <- varHaz$time*Year
    
    time.format <- "days" # Fix argument
  }
  else if(time.format == "months"){
    if(max(msfit.obj[[1]]$time) > 600) warning("Your max time in the data is more than 600 months. If time is not stored in months, please use argument time.format. \n")
    
    data$Tstart <- data$Tstart*Month
    data$Tstop <- data$Tstop*Month
    data$time <- data$time*Month
    
    Haz$time <- Haz$time*Month
    if(exists("varHaz")) varHaz$time <- varHaz$time*Month
    
    time.format <- "days" # Fix argument
  }
  else{
    stop("Argument time.format should take values in c('days', 'years', 'months').")
  }
  
  # Define new objects and transMat:
  Haz_new <- Haz[0,]
  if(length(msfit.obj) == 3){
    trans <- msfit.obj[[3]]
  }
  else{
    trans <- msfit.obj[[2]]
  }
  
  # Check split.transitions argument value:
  if(missing(split.transitions)) stop("Please define split.transitions.")
  else{
    if(class(split.transitions) == "numeric"){
      if(!all(split.transitions %in% 1:max(trans, na.rm=TRUE))) stop("Invalid transitions used inside argument split.transitions.")
    }
    else stop("Argument split.transitions expects values of class numeric")
    
    # Check intermediate states:
    for(i in split.transitions){
      tmp_state <- which(trans==i, arr.ind=TRUE)[2]
      if(!all(is.na(trans[tmp_state,]))) stop("You've listed an intermediate transition for which the hazard would have to be split in excess and population hazard. Population mortality tables for intermediate events haven't been implemented in this function. Please include only transitions that go to death states in the split.transitions argument.")
    }
  }
  
  # Get new transition matrix:
  trans_new <- modify_transMat(trans, split.transitions)
  
  # Get all transitions:
  transitions <- which(!is.na(trans), arr.ind = TRUE)
  link_trans <- list()
  
  ####################
  # 2. Prepare hazards
  ####################
  
  # We go through all original transitions and match them
  # to the ones in the new transMat (trans_new).
  # We then obtain the new hazards.
  
  for(i in 1:nrow(transitions)){
    # The transition in trans:
    trans_1 <- trans[transitions[i,1], transitions[i,2]]
    
    # We deal differently based on the type of transition
    # (whether we have to split the transition or not):
    if(!(trans_1 %in% split.transitions)){
      # The adequate transition in trans_new:
      trans_2 <- trans_new[rownames(trans)[transitions[i,1]],
                           colnames(trans)[transitions[i,2]]]
      # Save the linkage:
      link_trans[[ trans_1 ]] <- trans_2
      
      # Copy the hazards from the msfit object:
      df_tmp <- msfit.obj[[1]][msfit.obj[[1]]$trans == trans_1,]
      df_tmp$trans <- trans_2
      
      # Save:
      Haz_new <- rbind(Haz_new, df_tmp)
    }
    else{
      # The adequate transition in trans_new:
      trans_2 <- c(trans_new[rownames(trans)[transitions[i,1]],
                             paste0(colnames(trans)[transitions[i,2]], ".p")],
                   trans_new[rownames(trans)[transitions[i,1]],
                             paste0(colnames(trans)[transitions[i,2]], ".e")])
      # Save the linkage:
      link_trans[[ trans_1 ]] <- trans_2
      
      # Prepare pop. and excess hazard objects:
      df_p <- Haz[Haz$trans == trans_1, ]
      df_e <- Haz[Haz$trans == trans_1, ]
      
      # Take the subset we need:
      df_subset <- data[(data$from == transitions[i,1]) &
                          (data$to == transitions[i,2]),]
      
      # Find first time, when a jump happens:
      wh_jump <- which.max(df_e$Haz>0)
      is_jump <- any(df_e$Haz>0)
      if(!is_jump){
        if(!(link_trans_ind == TRUE & substitution == FALSE)){ # If it's not called when bootstrapping
          stop(paste0("There are no events occurring in transition ", trans_1, ". Please remove it from the split.transitions argument."))
        }
      }
      
      if(nrow(df_subset)==0 | nrow(df_p)==0){
        if(nrow(df_p)>0){
          df_p$trans <- trans_2[1]
          df_e$trans <- trans_2[2]
        }
      } else{
        # Calculate hazards:
        Hazs <- haz_function(Surv(Tstop, status)~1, data = df_subset, ratetable = ratetable,
                             rmap = rmap, 
                             add.times = df_p$time, include.all.times = FALSE)
        
        # Take care of late times, if present:
        # if(!all(df_p$time %in% Hazs$time)){
        #   additional_times <- sort(df_p$time[!(df_p$time %in% Hazs$time)])
        #   
        #   Hazs$time <- c(Hazs$time, additional_times)
        #   Hazs$haz.pop <- c(Hazs$haz.pop, rep(Hazs$haz.pop[length(Hazs$haz.pop)], length(additional_times)))
        #   Hazs$haz.excess <- c(Hazs$haz.excess, rep(Hazs$haz.excess[length(Hazs$haz.excess)], length(additional_times)))
        # }
        
        # Calculate hazards at the wanted times
        wh <- which(Hazs$time %in% (df_p$time))
        wh_l <- c(NA, wh[1:(length(wh)-1)])+1
        wh_l[1] <- 1
        
        haz.pop <- sapply(1:length(wh), 
                          function(x) sum(Hazs$haz.pop[wh_l[x]:wh[x]]))
        haz.excess <- sapply(1:length(wh),
                             function(x) sum(Hazs$haz.excess[wh_l[x]:wh[x]]))
        
        # Checks:
        # plot(1:length(haz.pop), (df_p$Haz[1:length(haz.pop)] - cumsum(haz.pop + haz.excess)), type="l")
        # plot(1:length(haz.pop), cumsum(haz.pop), type="l") 
        # plot(1:length(haz.pop), cumsum(haz.excess), type="l")
        
        # Population hazards:
        df_p$trans <- trans_2[1]
        df_p$Haz <- cumsum(haz.pop)
        
        # Excess hazards:
        df_e$trans <- trans_2[2]
        # df_e$Haz <- cumsum(haz.excess)
        df_e$Haz <- df_e$Haz - df_p$Haz # calculate like this (L_O - L_P), since it's more numerically stable.
        # Check: plot(1:length(haz.pop), (Haz[Haz$trans == trans_1, "Haz"] - df_e$Haz - df_p$Haz), type="l")
        
        # Make the times fully equal as in the msfit object:
        old_times <- msfit.obj[[1]]$time[msfit.obj[[1]]$trans == trans_1]
        df_p$time <- old_times
        df_e$time <- old_times  
        
        # Check if cum. excess haz<0 after first event time:
        if(wh_jump<nrow(df_e) & substitution & any(df_e$Haz[(wh_jump+1):nrow(df_e)]<0)){
          warning(paste0("The relative survival assumption that the observed hazard can be split in population and excess components might not hold for transition ", trans_1, ". Please check the estimated hazards. Consider removing transition ", trans_1, " from the split.transitions argument. \n"))
        }
      }
      
      # Save:
      Haz_new <- rbind(Haz_new, df_p, df_e)
    }
  }
  ordering <- order(Haz_new[,"trans"])
  Haz_new <- Haz_new[ordering,,drop=FALSE]
  rownames(Haz_new) <- 1:nrow(Haz_new)
  
  ###########################
  # 3. Prepare (co)variances:
  ###########################
  
  # Covariances, if needed:
  if(length(msfit.obj) == 3){
    
    if(sum(var.pop.haz %in% c('fixed', 'bootstrap', 'both'))==0){
      stop("Argument var.pop.haz should take value 'fixed', 'bootstrap' or 'both'.")
    }
    if(var.pop.haz == "fixed" | var.pop.haz == "both"){
      # Only adjust the original covariance matrix for the new transition matrix model:
      varHaz_new <- varHaz.fixed(varHaz, link_trans, varHaz_original=msfit.obj[[2]])
    }
    if(var.pop.haz == "bootstrap" | var.pop.haz == "both"){
      
      set.seed(seed)
      # Prepare object with all times:
      all_times <- Haz_new[, c("trans", "time")]
      colnames(all_times)[2] <- "Tstop"
      
      # Do the bootstrap:
      haz_boot <- msboot.relsurv(theta=msboot.relsurv.boot, data=data_orig, id=colnames(data_orig)[1],
                                 B=B, verbose=0, transmat=trans, all_times = all_times,
                                 split.transitions = split.transitions, rmap = rmap, time.format = time.format.orig, 
                                 boot_orig_msfit=FALSE, ratetable=ratetable, add.times=add.times
      )
      
      # Aggregate results:
      haz_boot_df <- do.call(rbind.data.frame, haz_boot)
      
      # Check if B is too small so that the variances are estimated
      # No. of bootstrap samples for every transition:
      no_samples <- table(unique(haz_boot_df[,c("trans", "b")])$trans)
      # If not enough samples for some transition or some transition not present in the end, stop:
      if(any(no_samples<2) | (!identical(as.integer(names(no_samples)), sort(as.integer(stats::na.omit(as.vector(trans_new))))))){
        stop("Not enough bootstrap samples so that the variances can be estimated. Please increase value in argument B.")
      }
      
      # CHECK:
      # pdf("boot_replicates.pdf")
      # ggplot(haz_boot_df %>% filter(trans %in% c(2,4)) %>% mutate(trans=ifelse(trans==2,"HCT->NRM.p","Rel->DaR.p")) %>% mutate(gr=paste0(trans, b), trans=factor(trans)), aes(time, Haz, group=gr, color=gr))+
      #   geom_line()+
      #   ylab("Hazards for different boot replicates")+
      #   facet_wrap(~trans, nrow = 1)+
      #   theme_bw()+
      #   theme(legend.position="none")
      # dev.off()
      
      haz_boot2 <-  aggregate(.~time+trans, data=subset(haz_boot_df, select=-b), FUN = stats::var)
      colnames(haz_boot2)[3] <- "varHaz"
      
      # browser()
      # Check bootstrap transitions:
      # juhu <- haz_boot_df %>% filter(trans==5, time<0.05)
      # juhu <- rbind(juhu %>% mutate(ju="b"),
      #               haz_boot2 %>% dplyr::select(time, Haz=varHaz, trans) %>% filter(trans==5, time<0.05) %>% mutate(b=1, ju="var"))
      # ggplot(juhu %>% mutate(b=factor(b)), aes(time, Haz, group=b, color=b))+
      #   geom_line(position=position_jitter(w=0, h=0.00))+
      #   facet_wrap(~ju, nrow=2, scales = "free_y")+
      #   theme_bw()+
      #   theme(legend.position = "none")
      
      # Add trans1, trans2:
      haz_boot3 <- haz_boot2
      haz_boot3$trans1 <- haz_boot3$trans
      colnames(haz_boot3)[2] <- "trans2"
      haz_boot3 <- haz_boot3[,c("time", "varHaz", "trans1", "trans2")]
      
      # Save objects:
      if(var.pop.haz == 'both') varHaz_new <- rbind(data.frame(varHaz_new, var.pop.haz="fixed"),
                                                    data.frame(haz_boot3, var.pop.haz="bootstrap"))
      else varHaz_new <- haz_boot3
      
      bootstrap_ind <- TRUE
      
      # CHECK
      # haz_boot <- haz_boot %>%
      #   group_by(time) %>%
      #   summarise(pstate3 = sd(pstate3, na.rm = TRUE),
      #             pstate4 = sd(pstate4, na.rm = TRUE),
      #             pstate5 = sd(pstate5, na.rm = TRUE),
      #             pstate6 = sd(pstate6, na.rm = TRUE)) %>%
      #   gather(state, prob, -time) %>%
      #   mutate(state = ifelse(state=="pstate3", "HCT->NRM.p",
      #                         ifelse(state=="pstate4", "HCT->NRM.e",
      #                                ifelse(state=="pstate5", "HCT->DaR.p", "HCT->DaR.e"))))
    }
    ordering <- order(varHaz_new[,c("trans1", "trans2")])
    varHaz_new <- varHaz_new[ordering,,drop=FALSE]
    varHaz_new <- stats::na.omit(varHaz_new)
    rownames(varHaz_new) <- 1:nrow(varHaz_new)
  }
  
  ####################
  # 4. Return objects:
  ####################
  
  # Save the new values:
  if(length(msfit.obj) == 3){
    res <- list(Haz=Haz_new,varHaz=varHaz_new,trans=trans_new)
  } else{
    res <- list(Haz=Haz_new,trans=trans_new)
  }
  
  # If link_trans and bootstrap replications have to be added:
  if(link_trans_ind) res$link_trans <- link_trans
  if(bootstrap_ind) res$Haz.boot <- haz_boot
  
  class(res) <- "msfit"
  return(res)
}
