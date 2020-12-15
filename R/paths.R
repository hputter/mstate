#' Find all possible trajectories through a given multi-state model
#' 
#' For a given multi-state model, specified through a transition matrix,
#' \code{paths} recursively finds all the possible trajectories or paths
#' through that multi-state starting from a specified state. DO NOT USE for
#' reversible or cyclic multi-state models.
#' 
#' The function is recursive. It starts in \code{start}, looks at what states
#' can be visited from \code{start}, and appends the results of the next call
#' to the current value (matrix). If the transition matrix contains loops, the
#' function will find infinitely many paths, so do not use \code{paths} for
#' reversible or cyclic multi-state models. A warning is not yet incorporated!
#' 
#' @param trans The transition matrix describing the multi-state model, see
#' \code{\link{msprep}}
#' @param start The starting state for the trajectories
#' @return A matrix, each row of which specifies a possible path through the
#' multi-state model.
#' @author Hein Putter <H.Putter@@lumc.nl>
#' @keywords array
#' @examples
#' 
#' tmat <- matrix(NA,5,5)
#' tmat[1,2:3] <- 1:2
#' tmat[1,5] <- 3
#' tmat[2,4:5] <- 4:5
#' tmat[3,4:5] <- 6:7
#' tmat[4,5] <- 8
#' paths(tmat)
#' paths(tmat, start=3)
#' 
#' @export paths
`paths` <- function(trans,start=1)
{
    ### Find all paths through a multi-state model from a starting state
    ### All direct transitions are specified through transition matrix
    ### trans
    ### Input:
    ###     start: numeric, specifying the starting state
    ###     trans: transition matrix
    ### Output:
    ###     matrix specifying all possible paths
    ### Details: recursive
    if (is.circular(trans)) stop("transition matrix is circular, so there will be infinitely many paths")
    nstates <- trans[start,]
    nstates <- which(!is.na(nstates))
    if (length(nstates)==0) ## i.e. in absorbing state
        return(matrix(start,1,1))
    else {
        fpmat <- startmat <- matrix(start,1,1)
        for (nstate in nstates) {
            fpmat <- my.rbind(fpmat, my.cbind2(start,Recall(trans,start=nstate)))
        }
    }
    return(fpmat)
}
