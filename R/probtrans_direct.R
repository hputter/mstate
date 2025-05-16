#' Given intensity matrices at each unique time, calculate transition probabilities
#' for a general multi-state model.
#' 
#' @description Given intensity matrices as output from \code{\link{get_intensity_matrices}},
#' calculate the transition probabilities. 
#'
#'
#' @param int_mat Output from \code{\link{get_intensity_matrices}}
#' @param predt Prediction time from which to make predictions
#' @param direction Should predictions be made from \code{predt} into the future ("fixed"),
#' or starting from time 0 until \code{predt} ("fixedhorizon").
#' @param as.df Should output be returned in \code{data.frame} format. Default is FALSE.
#' Setting this to TRUE returns exactly the same output as \code{\link[mstate:probtrans]{probtrans}}.
#'
#'
#' @details For more details, see \code{\link[mstate:probtrans]{probtrans}}.
#' 
#' @author Daniel Gomon
#' 
#' @return With \code{as.df = TRUE}, return exactly what \code{\link[mstate:probtrans]{probtrans}}
#' does. With \code{as.df = FALSE}, return an array instead, where the 3D dimension 
#' can be viewed as the list elements.
#'
#' @keywords internal
#' @noRd
#'
#'

probtrans_D <- function(int_mat, predt, direction = c("forward", "fixedhorizon"), as.df = FALSE){
  
  
  #Get only the matrices and times relevant for our probabilities
  if(direction == "forward"){
    relevant_idx <- predt < int_mat$unique_times
    relevant_matrices <- int_mat$intensity_matrices[, , relevant_idx, drop = FALSE ]
    times <- c(predt, int_mat$unique_times[relevant_idx])
  } else{
    relevant_idx <- predt >= int_mat$unique_times
    relevant_matrices <- int_mat$intensity_matrices[, , relevant_idx, drop = FALSE ]
    times <- int_mat$unique_times[relevant_idx]
    #DOES NOT CONTAIN predt if predt is not in unique_times!!
  }
  n_states <- dim(relevant_matrices)[1]
  n_matrices <- dim(relevant_matrices)[3]
  #We do not accept prediction times on boundaries or outside study time.
  if(n_matrices == 0){
    stop("Please specify prediction time 'predt' within study time (boundaries also not allowed).
         Transition probabilities could not be calculated.")
  }
  
  #We are certain to stay in initial state at initial time point ("forward")
  #and final time point ("fixedhorizon")
  P <- diag(n_states)
  
  
  #Forward direction
  if(direction == "forward"){
    #Extra column to keep track of times
    out <- array(0, dim = c(n_states, n_states + 1, n_matrices + 1))
    #For the first transition time, we know we must stay in the same state at predt coming from predt
    for(j in 1:n_states){
      out[j, 1, ] <- times
    }
    out[, 2:(n_states +1) , 1] <- P
    #For other transition times, multiply the corresponding matrices
    for(i in 1:n_matrices){
      P <- P %*% relevant_matrices[, , i ]
      out[, 2:(n_states +1) , i+1] <- P
    }  
  } else{ #Backward direction
    #If prediction time is not in unique_times, we need to create extra entry in probtrans
    if(predt %in% times){
      out <- array(0, dim = c(n_states, n_states + 1, n_matrices + 1))
      for(j in 1:n_states){
        out[j, 1, ] <- c(0, times)
      }
      out[, 2:(n_states + 1), n_matrices + 1] <- P
    } else{
      out <- array(0, dim = c(n_states, n_states + 1, n_matrices + 2))
      for(j in 1:n_states){
        out[j, 1, ] <- c(0, times, predt)
      }
      #Sure to stay in state at final time points
      out[, 2:(n_states + 1), n_matrices + 1] <- P
      out[, 2:(n_states + 1), n_matrices + 2] <- P
    }
    
    #For other transition times, multiply the corresponding matrices
    for(i in 1:n_matrices){
      P <- relevant_matrices[, , n_matrices + 1 - i ] %*% P
      out[, 2:(n_states +1) , n_matrices + 1 - i ] <- P
    }
  }
  if(as.df){ #get output exactly as mstate::probtrans()
    #Store in list of data frames
    res <- vector(mode = "list", length = n_states)
    for(s in 1:n_states){
      res[[s]] <- as.data.frame(t(out[s, , ]))
      colnames(res[[s]]) <- c("time", paste("pstate", 1:n_states, sep = ""))
    }  
  } else{ #Output 3D array instead
    #Suppose a <- mstate::probtrans() output
    #Then a[[1]] == res[, , 1]
    res <- aperm(out, perm = c(3, 2 , 1))
  }
  
  return(res)
}


#' Given intensity matrices at each unique time, calculate transition probabilities
#' for a general multi-state model.
#' 
#' @description Given intensity matrices as output from \code{\link{get_intensity_matrices}},
#' calculate the transition probabilities, allowing for a cut-off time to be specified. Basically 
#' a clone of \code{\link{probtrans_D}} but allows to specify cutoff.
#'
#'
#' @param int_mat Output from \code{\link{get_intensity_matrices}}
#' @param predt Prediction time from which to make predictions
#' @param cutoff Cutoff for prediction times.
#' @param direction Should predictions be made from \code{predt} into the future ("fixed") until \code{cutoff},
#' or starting from \code{cutoff} until \code{predt} ("fixedhorizon").
#' @param as.df Should output be returned in \code{data.frame} format. Default is FALSE.
#' Setting this to TRUE returns exactly the same output as \code{\link[mstate:probtrans]{probtrans}}.
#'
#'
#' @details For more details, see \code{\link[mstate:probtrans]{probtrans}}.
#' 
#' @author Daniel Gomon
#' 
#' @return With \code{as.df = TRUE}, return exactly what \code{\link[mstate:probtrans]{probtrans}}
#' does. With \code{as.df = FALSE}, return an array instead, where the 3D dimension 
#' can be viewed as the list elements.
#'
#' @keywords internal
#' @noRd
#'
#'

probtrans_C <- function(int_mat, predt, cutoff, direction = c("forward", "fixedhorizon"), as.df = FALSE){
  
  #Only difference between probtrans_C and probtrans_D is that
  #probtrans_C has the cutoff argument, allowing to select smaller time windows
  #Cutoff must be larger/smaller than predt when direction is forward/fixedhorizon
  
  #Get only the matrices and times relevant for our probabilities
  if(direction == "forward"){
    relevant_idx <- predt < int_mat$unique_times & int_mat$unique_times <= cutoff
    relevant_matrices <- int_mat$intensity_matrices[, , relevant_idx, drop = FALSE ]
    times <- c(predt, int_mat$unique_times[relevant_idx])
  } else{
    relevant_idx <- cutoff < int_mat$unique_times & int_mat$unique_times <= predt
    relevant_matrices <- int_mat$intensity_matrices[, , relevant_idx, drop = FALSE ]
    times <- int_mat$unique_times[relevant_idx]
    #DOES NOT CONTAIN predt if predt is not in unique_times!!
  }
  n_states <- dim(relevant_matrices)[1]
  n_matrices <- dim(relevant_matrices)[3]
  #We do not accept prediction times on boundaries or outside study time.
  if(n_matrices == 0){
    stop("Please specify prediction time 'predt' within study time (boundaries also not allowed).
         Transition probabilities could not be calculated.")
  }
  
  #We are certain to stay in initial state at initial time point ("forward")
  #and final time point ("fixedhorizon")
  P <- diag(n_states)
  
  
  #Forward direction
  if(direction == "forward"){
    #Extra column to keep track of times
    out <- array(0, dim = c(n_states, n_states + 1, n_matrices + 1))
    #For the first transition time, we know we must stay in the same state at predt coming from predt
    for(j in 1:n_states){
      out[j, 1, ] <- times
    }
    out[, 2:(n_states +1) , 1] <- P
    #For other transition times, multiply the corresponding matrices
    for(i in 1:n_matrices){
      P <- P %*% relevant_matrices[, , i ]
      out[, 2:(n_states +1) , i+1] <- P
    }  
  } else{ #Backward direction "fixedhorizon"
    #If prediction time is not in unique_times, we need to create extra entry in probtrans
    if(predt %in% times){
      out <- array(0, dim = c(n_states, n_states + 1, n_matrices + 1))
      for(j in 1:n_states){
        out[j, 1, ] <- c(cutoff, times)
      }
      out[, 2:(n_states + 1), n_matrices + 1] <- P
    } else{
      out <- array(0, dim = c(n_states, n_states + 1, n_matrices + 2))
      for(j in 1:n_states){
        out[j, 1, ] <- c(cutoff, times, predt)
      }
      #Sure to stay in state at final time points
      out[, 2:(n_states + 1), n_matrices + 1] <- P
      out[, 2:(n_states + 1), n_matrices + 2] <- P
    }
    
    #For other transition times, multiply the corresponding matrices
    for(i in 1:n_matrices){
      P <- relevant_matrices[, , n_matrices + 1 - i ] %*% P
      out[, 2:(n_states +1) , n_matrices + 1 - i ] <- P
    }
  }
  if(as.df){ #get output exactly as mstate::probtrans()
    #Store in list of data frames
    res <- vector(mode = "list", length = n_states)
    for(s in 1:n_states){
      res[[s]] <- as.data.frame(t(out[s, , ]))
      colnames(res[[s]]) <- c("time", paste("pstate", 1:n_states, sep = ""))
    }  
  } else{ #Output 3D array instead
    #Suppose a <- mstate::probtrans() output
    #Then a[[1]] == res[, , 1]
    res <- aperm(out, perm = c(3, 2 , 1))
  }
  
  return(res)
}








#' Given an \code{msfit} object, extract intensity matrices in 3D array form
#' 
#' @description Transform an \code{\link[mstate:msfit]{msfit}} object into 
#' a 3D array representing the intensity matrix over time. The dimensions represent:
#' \describe{
#'   \item{\code{1st}:}{States from which a transition takes place;}
#'   \item{\code{2nd}:}{States to which a transition takes place;}
#'   \item{\code{3rd}:}{Unique times in the \code{\link[mstate:msfit]{msfit}} 
#'   object;}
#' } Note that for an intensity matrix the rows must sum to 1, so fixing the 
#' 3rd dimension the diagonal elements are one minus the row sums.
#' 
#' @param object An \code{\link[mstate:msfit]{msfit}} object.
#' 
#' @author Daniel Gomon
#' 
#' @keywords internal
#' @noRd
#' 
#' 

get_intensity_matrices <- function(object){
  
  tmat2 <- mstate::to.trans2(object$trans)
  
  #object$Haz must be sorted according to trans, then time!!
  unique_times <- unique(object$Haz$time)
  n_times <- length(unique_times)
  n_states <- nrow(object$trans)
  n_trans <- nrow(tmat2)
  
  #Initialize output matrices
  intensity_matrices <- array(0, dim = c(n_states, n_states, n_times))
  
  #Determine intensities
  for(k in 1:n_trans){
    t_obj <- object$Haz[object$Haz$trans == k,]
    intensities <- diff(c(0, t_obj$Haz))
    from <- tmat2$from[k]
    to <- tmat2$to[k]
    
    #Assign determined intensities to corresponding matrix
    intensity_matrices[from, to, ] <- intensities
  }
  
  for(i in 1:n_times){
    diag(intensity_matrices[, , i]) <- 1 - rowSums(intensity_matrices[, , i])
  }
  out <- list(intensity_matrices = intensity_matrices,
              unique_times = unique_times)
  return(out)
}






