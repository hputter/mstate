#' @export
`to.trans2` <- function(trans)
{
  dm <- dim(trans)
  if (dm[1] != dm[2]) stop("transition matrix should be square")
  S <- dm[1]
  mx <- max(trans, na.rm=TRUE)
  res <- matrix(NA, mx, 3)
  res[, 1] <- 1:mx
  transvec <- as.vector(trans)
  for (i in 1:mx) {
    idx <- which(transvec==i)
    res[i, 2:3] <- c((idx-1) %% S + 1, (idx-1) %/% S + 1)
  }
  res <- data.frame(res)
  names(res) <- c("transno", "from", "to")
  res$from[res$from==0] <- S
  statesfrom <- dimnames(trans)[[1]]
  if (is.null(statesfrom)) statesfrom <- 1:S
  statesto <- dimnames(trans)[[2]]
  if (is.null(statesto)) statesto <- 1:S
  res$fromname <- statesfrom[res$from]
  res$toname <- statesto[res$to]
  res$transname <- paste(res$fromname, res$toname, sep=" -> ")
  return(res)
}

#' @export
`trans2Q` <- function(trans)
{
  K <- nrow(trans)
  P <- trans
  P[!is.na(P)] <- 1
  P[is.na(P)] <- 0
  diag(P) <- 1
  k <- 1
  # deb(k, method="cat")
  Pk <- P
  diag(Pk) <- 0
  # deb(Pk)
  Pkprev <- Pk
  Q <- Pk
  for (k in 2:K) {
    # deb(k, method="cat")
    Pk <- Pk %*% P
    Pk[Pk > 1] <- 1
    # deb(Pk)
    # deb(Pk - Pkprev)
    Q <- Q + k * (Pk - Pkprev)
    # deb(Q)
    Pkprev <- Pk
  }
  Q
}

#' @export
`absorbing` <- function(trans)
{
  Q <- trans2Q(trans)
  wh <- which(apply(Q, 1, sum) == 0)
  wh
}

#' @export
`is.circular` <- function(trans)
{
  Q <- trans2Q(trans)
  return(any(diag(Q)>0))
}
