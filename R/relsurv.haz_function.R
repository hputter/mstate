#' Helper function that calculates excess and population hazards for a given transition
#' 
#' A function that calculates the excess and population hazards
#' for a given transition. Code is based on function rs.surv from
#' the relsurv package. 
#' @param formula A non-parametric Surv-based formula, e.g. Surv(times, status)~1.
#' @param data A subset of the msprep object (dataset) where there's 
#' only data for the chosen transition.
#' @param ratetable A table of event rates, organized as a ratetable object, such as slopop.
#' @param na.action A missing-data filter function, applied to the model.frame, after any subset argument has been used. Default is options()$na.action.
#' @param add.times Additional times at which the hazards should be evaluated.
#' @param rmap An optional list to be used if the variables are not organized and named in the same way as in the ratetable object
#' @param include.all.times Should hazards be evaluated at all times in seq(minimum time, maximum time, by=1). Default is FALSE

#' @seealso \code{\link{msfit.relsurv}}
#' @export 
`haz_function` <- function(formula = formula(data), data, ratetable = relsurv::slopop,
                         na.action, add.times, rmap, include.all.times=FALSE
){
  
  call <- match.call()
  if(!missing(rmap))  rmap <- as.call(rmap)
  
  # if(data$trans[1]==3){ browser()}
  
  # Adjust the dataset:
  # status_indicator <- as.character(terms(formula)[[2]][3]) # find the name of the status indicator
  # data[data[, status_indicator] != failcode, status_indicator] <- 0 # censor those that don't go to the desired state
  # data[data[, status_indicator] == failcode, status_indicator] <- 1 # make sure the status indicator is equal to 1
  
  data_yi <- data # Save original data
  
  # Prepare rform:
  if(nrow(data)==0) browser()
  rform <- relsurv:::rformulate(formula,data,ratetable,na.action,rmap)
  data <- rform$data # dataset
  
  # Adjust year/age for left truncation:
  # if(left.truncation){
  #   
  #   Tstart <- round(data_yi$Tstart)
  #   # data:
  #   rform$data$year <- rform$data$year + Tstart
  #   rform$data$age <- rform$data$age + Tstart
  #   
  #   # R:
  #   rform$R[,"year"] <- rform$R[,"year"] + Tstart
  #   rform$R[,"age"] <- rform$R[,"age"] + Tstart
  # }
  
  # Check covariates:
  p <- rform$m
  if (p > 0) 
    stop("There shouldn't be any covariates in the formula. This function gives non-parametric estimates of the hazards.")
  else data$Xs <- rep(1, nrow(data)) #if no covariates, just put 1
  
  out <- NULL
  out$n <- table(data$Xs) #table of strata
  out$time <- out$n.risk <- out$n.event <- out$n.censor <- haz.excess <- haz.pop <- out$std.err <- NULL
  
  kt <- 1 # the only stratum
  inx <- which(data$Xs == names(out$n)[kt]) #individuals within this stratum
  tis <- sort(unique(rform$Y[inx])) #unique times
  
  if(include.all.times){
    # Include all times between 1 and maximum time (by 1 day):
    all.steps <- seq(min(rform$Y[inx]), max(rform$Y[inx]), by = 1)
    tis <- sort(union(rform$Y[inx],as.numeric(all.steps)))
  }
  
  # Include add.times, if needed:
  if(!missing(add.times)){
    # add.times <- pmin(as.numeric(add.times),max(rform$Y[inx]))
    tis <- sort(union(tis,as.numeric(add.times)))	#1-day long intervals used - to take into the account the continuity of the pop. part
  }
  
  # Left-truncation in R:
  # if(left.truncation){
  #   # Prepare at-risk matrix:
  #   mat <- lapply(1:nrow(data_yi), function(x) ifelse((data_yi$Tstart[x] < tis) & (tis <= data_yi$Tstop[x]), 1, 0))
  #   mat2 <- matrix(unlist(mat), nrow = nrow(data_yi), byrow = TRUE)
  #   # The sum of the individual at-risk processes:
  #   yi_left <- colSums(mat2)
  #   yi_left[yi_left == 0] <- Inf
  #   
  #   # Calculate yidli (part of pop. hazard) for left-truncated data.
  #   # Prepare matrix:
  #   max_year <- as.character(max(as.integer(colnames(ratetable)))) # Maximum year in ratetable
  #   data$year_helper <- as.Date(data$year, origin="1960-01-01") # Helper column
  #   
  #   
  #   mat4 <- lapply(1:nrow(data_yi), function(i) lt_fun(i, mat2, tis, ratetable, data, max_year))
  #   mat5 <- matrix(unlist(mat4), nrow = nrow(data), byrow = TRUE)
  #   
  #   # mat5 <- lt_fun2(mat2, tis, ratetable, data, max_year)
  #   
  #   # yidli_by_hand <- colSums(mat3) # do the summation across individuals
  #   yidli_by_hand <- colSums(mat5) # do the summation across individuals
  # }
  # data$year_helper <- NULL # remove helper column
  
  # Calculate the values for each interval of time:
  temp <- relsurv:::exp.prep(rform$R[inx,,drop=FALSE],rform$Y[inx],rform$ratetable,rform$status[inx],times=tis,fast=TRUE, cmp=FALSE, ys=data_yi$Tstart)	
  
  # Fix at-risk process, if needed:
  temp$yi[temp$yi==0] <- Inf
  # for(ii in 1:length(temp$yi)){
  #   if(temp$yi[ii]==0){
  #     temp$yi[ii] <- Inf
  #   }
  #   else{
  #     break
  #   }
  # }
  
  out$time <- c(out$time, tis)						#add times
  out$n.risk <- c(out$n.risk, temp$yi)					#add number at risk for each time
  out$n.event <- c(out$n.event, temp$dni)					#add number of events for each time
  out$n.censor <- c(out$n.censor,  c(-diff(temp$yi),temp$yi[length(temp$yi)]) - temp$dni) 	#add number of censored for each time
  
  # Calculate hazards:
  haz.excess <- temp$dni/temp$yi - temp$yidli/temp$yi
  haz.pop <- temp$yidli/temp$yi
  
  out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))  #standard error on each interval
  
  out$haz.excess <- c(out$haz.excess,haz.excess)
  out$haz.pop <- c(out$haz.pop,haz.pop)
  
  out$n <- as.vector(out$n)
  out$call <- call
  out
}