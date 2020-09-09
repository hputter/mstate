#' Define transition matrix for multi-state model
#' 
#' Define transition matrices for multi-state model. Specific functions for
#' defining such transition matrices are pre-defined for common multi-state
#' models like the competing risks model and the illness-death model.
#' 
#' If \code{names} is missing, the names \code{"eventfree"}, \code{"cause1"},
#' etcetera are assigned in \code{trans.comprisk}, or \code{"healthy"},
#' \code{"illness"}, \code{"death"} in \code{trans.illdeath}.
#' 
#' @aliases transMat trans.illdeath trans.comprisk
#' @param x List of possible transitions; x[[i]] consists of a vector of state
#' numbers reachable from state i
#' @param K The number of competing risks
#' @param names A character vector containing the names of either the competing
#' risks or the states in the multi-state model specified by the competing
#' risks or illness-death model. \code{names} should have the same length as
#' the list \code{x} (for \code{transMat}), or either \code{K} or \code{K}+1
#' (for \code{trans.comprisk}), or 3 (for \code{trans.illdeath})
#' @return A transition matrix describing the states and transitions in the
#' multi-state model.
#' @author Steven McKinney <smckinney@@bccrc.ca>; Hein Putter
#' <H.Putter@@lumc.nl>
#' @keywords array
#' @examples
#' 
#' transMat(list(c(2, 3), c(), c(1, 2)),
#' 	names = c("Disease-free", "Death", "Relapsed"))
#' tmat <- transMat(x = list( c(2, 3), c(1, 3), c() ),
#'                  names = c("Normal", "Low", "Death"))
#' tmat
#' transListn <- list("Normal" = c(2, 3), "Low" = c(1, 3), "Death" = c())
#' transMat(transListn)
#' trans.comprisk(3)
#' trans.comprisk(3,c("c1","c2","c3"))
#' trans.comprisk(3,c("nothing","c1","c2","c3"))
#' trans.illdeath()
#' trans.illdeath(c("nothing","ill","death"))
#' 
#' @export transMat
transMat <- function(x, names) {
  ## transMat:  produce transition matrix for use in package 'mstate'.
  ## Arguments:
  ## x:  List of possible transitions.
  ##     x[[i]] consists of a vector of state numbers
  ##     reachable from state i.
  ## names: Character vector of state names, having the same length
  ##        as x.
  ## Example:  States 1 and 2 are reachable from one
  ##     another.  State 3 is absorbing, reachable from
  ##     states 1 and 2.
  ##     transMat( x = list( c(2, 3), c(1, 3), c() ) )
  if ( !is.list(x) ) stop("x must be a list")
  ns <- length(x) ## number of states
  tmat <- matrix(NA, nrow = ns, ncol = ns) ## transition matrix
  if ( missing(names) ) {
    if ( !is.null( base::names(x) ) ) {
      namesList <- list(from = base::names(x), to = base::names(x))
    } else {
      namesList <- list(from = paste("State", seq(nrow(tmat))),
                        to = paste("State", seq(nrow(tmat))))
    }
  } else {
    if ( length(names) != ns ) stop("length of 'names' must equal length of 'x'")
    namesList <- list(from = names, to = names)
  }
  idxmat <- cbind(unlist(lapply(seq(ns),
                                function(i, y){
                                  rep(i, length(y[[i]]))}, x)),
                  unlist(x))
  if ( max(idxmat) > ns )
    stop("Largest state in transition list exceeds number of states")
  tmat[idxmat] <- seq(nrow(idxmat))
  dimnames(tmat) <- namesList
  tmat
}
