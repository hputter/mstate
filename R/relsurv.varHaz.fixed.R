#' Upgrade the varHaz object
#'
#' A function that upgrades varHaz from the msfit object where the variances are
#' estimated using the Greenwood estimator; it is further assumed that variances for the population
#' hazards are equal to zero.
#' @param varHaz The varHaz object (present in a msfit object).
#' @param link_trans A list that gives the linkage between the original and upgraded 
#' transition matrix.
#' @param varHaz_original The original varHaz object from msfit (without the eventual time conversion).
#' @return Return the upgraded varHaz object containing variances for the split transitions.
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' 
#' @export
`varHaz.fixed` <- function(varHaz, link_trans, varHaz_original){
  
  # Define new varHaz object:
  varHaz_new <- varHaz[0,]
  
  # Find all covariance combinations:
  var_trans <- unique(varHaz[,c("trans1","trans2")])
  
  # For every combination write down the (co)variance:
  for(i in 1:nrow(var_trans)){
    
    trans1 <- var_trans$trans1[i]
    trans2 <- var_trans$trans2[i]
    
    # Find the adequate transitions in trans_new:
    trans1_new <- link_trans[[ trans1 ]]
    trans2_new <- link_trans[[ trans2 ]]
    
    varHaz_tmp <- varHaz[(varHaz$trans1 == trans1) &
                           (varHaz$trans2 == trans2), ]
    
    if((length(trans1_new) == 1) & (length(trans2_new) == 1)){
      # Replace suitable tranistions and add matrix:
      varHaz_tmp$trans1 <- trans1_new
      varHaz_tmp$trans2 <- trans2_new
      varHaz_tmp$time <- varHaz_original$time[(varHaz_original$trans1 == trans1) &
                                                (varHaz_original$trans2 == trans2)]
      varHaz_new <- rbind(varHaz_new, varHaz_tmp)
    }
    else if(length(trans1_new) == 1){
      # Replace suitable tranistions and add matrices:
      varHaz_tmp$trans1 <- trans1_new
      varHaz_tmp$trans2 <- trans2_new[2]
      varHaz_tmp$time <- varHaz_original$time[(varHaz_original$trans1 == trans1) &
                                                (varHaz_original$trans2 == trans2)]
      varHaz_tmp2 <- varHaz_tmp
      
      varHaz_tmp$trans2 <- trans2_new[1]
      varHaz_tmp$varHaz <- 0
      varHaz_new <- rbind(varHaz_new, varHaz_tmp, varHaz_tmp2)
    }
    else if(length(trans2_new) == 1){
      # Replace suitable tranistions and add matrices:
      varHaz_tmp$trans1 <- trans1_new[2]
      varHaz_tmp$trans2 <- trans2_new
      varHaz_tmp$time <- varHaz_original$time[(varHaz_original$trans1 == trans1) &
                                                (varHaz_original$trans2 == trans2)]
      varHaz_tmp2 <- varHaz_tmp
      
      varHaz_tmp$trans1 <- trans1_new[1]
      varHaz_tmp$varHaz <- 0
      varHaz_new <- rbind(varHaz_new, varHaz_tmp, varHaz_tmp2)
    }
    else{
      # Replace suitable tranistions and add matrices:
      varHaz_tmp$trans1 <- trans1_new[2]
      varHaz_tmp$trans2 <- trans2_new[2]
      varHaz_tmp$time <- varHaz_original$time[(varHaz_original$trans1 == trans1) &
                                                (varHaz_original$trans2 == trans2)]
      varHaz_tmp2 <- varHaz_tmp
      
      varHaz_tmp$trans1 <- trans1_new[1]
      varHaz_tmp$trans2 <- trans2_new[1]
      varHaz_tmp$varHaz <- 0
      varHaz_new <- rbind(varHaz_new, varHaz_tmp)
      
      varHaz_tmp$trans1 <- trans1_new[1]
      varHaz_tmp$trans2 <- trans2_new[2]
      varHaz_new <- rbind(varHaz_new, varHaz_tmp)
      
      if(!identical(trans1_new, trans2_new)){
        varHaz_tmp$trans1 <- trans1_new[2]
        varHaz_tmp$trans2 <- trans2_new[1]
        varHaz_new <- rbind(varHaz_new, varHaz_tmp)
      }
      varHaz_new <- rbind(varHaz_new, varHaz_tmp2)
    }
  }
  
  return(varHaz_new)
}