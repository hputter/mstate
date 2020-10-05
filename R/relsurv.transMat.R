#' Upgrade the transMat object for the multi-state/relsurv setting.
#' 
#' A function that upgrades the transMat object so that the population
#' and excess-related transitions are included in the transition matrix.
#' @param trans The original transition matrix (usually generated using function
#' transMat from mstate). Also often present in the msfit object.
#' @param split.transitions The transitions that should be split.
#' @return An upgraded transition matrix that contains the population and
#' excess transitions.
#' @author Damjan Manevski \email{damjan.manevski@@mf.uni-lj.si}
#' @seealso \code{\link{transMat}}
#' @examples
#' 
#' # transition matrix for illness-death model
#' trans <- transMat(list(c(2,3),c(4), c(), c()), 
#'       names = c("Alive", "Relapse","Non-relapse mortality", "Death after relapse"))
#' split.transitions <- c(2,3)
#' modify_transMat(trans, split.transitions) 
#' 
#' @export 
`modify_transMat` <- function(trans, split.transitions){
  
  # Define objects:
  trans_tmp <- trans
  companions <- c()
  split.transitions <- sort(split.transitions)
  split.states <- sapply(split.transitions, function(x) which(trans == x, arr.ind=TRUE)[,2])
  
  # The function works in three steps:
  
  # 1. Add pop. and excess states and the additional transitions in trans:
  for(j in 1:ncol(trans)){
    
    check_trans <- trans[!is.na(trans[,j]),j] %in% split.transitions
    
    # If state j shouldn't be split, or all transitions going to j 
    # have to be split, then:
    if(sum(check_trans) %in% c(0, length(check_trans))){
      
      add_col <- TRUE
      for(i in 1:nrow(trans)){
        if(!is.na(trans[i,j])){
          if(j %in% split.states){
            if(add_col){
              # Add column:
              trans_tmp <- cbind(trans_tmp, NA)
              colnames(trans_tmp)[ncol(trans_tmp)] <- paste0(colnames(trans_tmp)[j], ".e")
              colnames(trans_tmp)[j] <- paste0(colnames(trans_tmp)[j], ".p")
              # Add row:
              trans_tmp <- rbind(trans_tmp, NA)
              rownames(trans_tmp)[nrow(trans_tmp)] <- paste0(rownames(trans_tmp)[j], ".e")
              rownames(trans_tmp)[j] <- paste0(rownames(trans_tmp)[j], ".p")
              # Save companions (.p and .e)
              companions[j] <- ncol(trans_tmp)
              add_col <- FALSE
            } 
            if(trans_tmp[i,j] %in% split.transitions){
              trans_tmp[i, ncol(trans_tmp)] <- trans_tmp[i,j] + 1
            }
          }
        }
      }
    }
    # Else if some transitions going to state j have to be split, others not:
    else{
      # First make new state where transitions are not splitted:
      # Add column:
      trans_tmp <- cbind(trans_tmp, trans_tmp[,j])
      last_col <- ncol(trans_tmp)
      colnames(trans_tmp)[last_col] <- paste0(colnames(trans_tmp)[j], "")
      # Which values have to be removed:
      remove_values <- trans_tmp[,last_col] %in% split.transitions
      trans_tmp[remove_values, last_col] <- NA
      # Add row:
      trans_tmp <- rbind(trans_tmp, NA)
      rownames(trans_tmp)[nrow(trans_tmp)] <- paste0(rownames(trans_tmp)[j], "")
      ###########################
      # Now add pop and excess states:
      
      # Which values have to be removed:
      remove_values <- !(trans_tmp[,j] %in% split.transitions)
      remove_values <- remove_values & !is.na(trans_tmp[,j]) # remove NAs
      trans_tmp[remove_values, j] <- NA
      # Add column:
      trans_tmp <- cbind(trans_tmp, NA)
      colnames(trans_tmp)[ncol(trans_tmp)] <- paste0(colnames(trans_tmp)[j], ".e")
      colnames(trans_tmp)[j] <- paste0(colnames(trans_tmp)[j], ".p")
      # Add row:
      trans_tmp <- rbind(trans_tmp, NA)
      rownames(trans_tmp)[nrow(trans_tmp)] <- paste0(rownames(trans_tmp)[j], ".e")
      rownames(trans_tmp)[j] <- paste0(rownames(trans_tmp)[j], ".p")
      # Save companions (.p and .e)
      companions[j] <- ncol(trans_tmp)
      trans_tmp[, ncol(trans_tmp)] <- trans_tmp[,j] + 1
    }
  }
  
  # 2. Order the columns/rows in a sensible way
  # (splitted states .p, .e. are put next to each other):
  comp_vek <- which(!is.na(companions))
  ordered_cols <- setdiff(1:ncol(trans_tmp), c(comp_vek, companions[comp_vek]))
  ordered_cols <- c(ordered_cols, unlist(lapply(comp_vek, function(x) c(x, companions[x]))))
  trans_tmp <- trans_tmp[ordered_cols, ordered_cols]
  
  # 3. Make the numbering of transitions uniform:
  val <- 1
  for(i in 1:nrow(trans_tmp)){
    for(j in 1:ncol(trans_tmp)){
      if(!is.na(trans_tmp[i,j])){
        trans_tmp[i,j] <- val
        val <- val+1
      }
    }
  }
  # Add from/to as dimension names (this gets lost when cbind/rbind is used):
  names(dimnames(trans_tmp)) <- names(dimnames(trans))
  
  return(trans_tmp)
}
