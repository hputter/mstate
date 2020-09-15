#' msdata to etm format
#' 
#' @param msdata Multi-state data in \code{msdata} format, as used in
#' \code{mstate}
#' @inheritParams etm2msdata 
#' 
#' @export
msdata2etm <- function(msdata, id, covs)
{
  if (missing(id)) id <- "id"
  msdata$tostat <- msdata$to * msdata$status
  aggr <- aggregate(msdata[, c("status", "tostat")],
                    by=list(id=msdata[, id], from=msdata$from, Tstart=msdata$Tstart, Tstop=msdata$Tstop),
                    FUN=max)
  names(aggr)[names(aggr)=="tostat"] <- "to"
  ord <- c(1, 3, 4, 2, 6)
  etm <- aggr[, ord]
  names(etm)[2:5] <- c("entry", "exit", "from", "to")
  etm$to[etm$to==0] <- 99
  etm <- etm[order(etm[, id], etm$entry), ]
  
  # Add covariates
  if (missing(covs)) covs <- NULL
  ncovs <- length(covs)
  if (ncovs > 0) {
    msdatacovs <- as.data.frame(msdata[, c("id", covs)])
    msdatacovs <- msdatacovs[!duplicated(msdatacovs$id), ]
    etm <- merge(etm, msdatacovs, by="id")
  }
  etm
}

#' Convert transition matrix from mstate to etm format
#' 
#' @param trans Transition matrix in \code{mstate} format
#' 
#' @export
trans2tra <- function(trans)
  return(!(is.na(trans)))


#' Converts between etm and msdata format
#' 
#' Converts multi-state data back and forth between etm and msdata formats.
#' Covariates have to be dealt with separately.
#' 
#' \code{msdata2etm} will convert from \code{msdata} format to \code{etm}
#' format; \code{etm2msdata} will convert from \code{etm} format to
#' \code{msdata} format. Both \code{msdata2etm} and \code{etm2msdata} work with
#' basic time-fixed covariates. Time-dependent covariates are not supported.
#' The function \code{msdata2etm} will work for transition-specific covariates,
#' but the result does not really make much sense when used in etm.
#' 
#' @aliases etm2msdata tra2trans
#' @param id Column name identifying the subject id
#' @param etmdata Multi-state data in \code{etm} format
#' @param tra Transition matrix in \code{etm} format
#' @param covs Vector of column names containing covariates to be included
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @keywords datagen
#' @examples
#' 
#' # Transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # Data in wide format, for transition 1 this is dataset E1 of
#' # Therneau & Grambsch (T&G)
#' tg <- data.frame(id=1:6,illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'                  dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'                  x1=c(1,1,1,0,0,0),x2=c(6:1))
#' # Data in long format using msprep
#' tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
#'                  data=tg,keep=c("x1","x2"),trans=tmat, id="id")
#' # Same thing in etm format
#' tra <- trans2tra(tmat)
#' tgetm <- msdata2etm(tglong, id="id")
#' tgetm <- msdata2etm(tglong, id="id", covs=c("x1", "x2")) # with covariates
#' # And back
#' etm2msdata(tgetm, id="id", tra=tra)
#' etm2msdata(tgetm, id="id", tra=tra, covs=c("x1", "x2")) # with covariates
#' @export
etm2msdata <- function(etmdata, id, tra, covs)
{
  nout <- apply(tra, 1, sum)
  trans <- tra2trans(tra)
  to <- apply(trans, 1, function(x) which(!is.na(x)))
  trans <- apply(trans, 1, function(x) x[!is.na(x)])
  if (missing(id)) id <- "id"
  msdata <- data.frame(id = rep(etmdata[, id], nout[etmdata$from]),
                       from = rep(etmdata$from, nout[etmdata$from]),
                       to = unlist(lapply(etmdata$from, function(x) to[[x]])),
                       trans = unlist(lapply(etmdata$from, function(x) trans[[x]])),
                       # to = unlist(lapply(etmdata$from, function(x) to[, x])),
                       # trans = unlist(lapply(etmdata$from, function(x) trans[, x])),
                       Tstart = rep(etmdata$entry, nout[etmdata$from]),
                       Tstop = rep(etmdata$exit, nout[etmdata$from]),
                       time = NA,
                       status = rep(etmdata$to, nout[etmdata$from]))
  msdata$time <- msdata$Tstop - msdata$Tstart
  msdata$status[msdata$status==msdata$to] <- -1 # these are going to be 1
  msdata$status[msdata$status>0] <- 0
  msdata$status <- -msdata$status

  # Add covariates
  if (missing(covs)) covs <- NULL
  ncovs <- length(covs)
  if (ncovs > 0) {
    etmdatacovs <- as.data.frame(etmdata[, c("id", covs)])
    etmdatacovs <- etmdatacovs[!duplicated(etmdatacovs$id), ]
    msdata <- merge(msdata, etmdatacovs, by="id")
  }
  
  class(msdata) <- c("msdata", "data.frame")
  attr(msdata, "trans") <- trans
  return(msdata) 
}

#' @export
tra2trans <- function(tra)
{
  ttrans <- t(matrix(as.numeric(tra), nrow(tra), ncol(tra)))
  ntr <- sum(tra)
  ttrans[ttrans==1] <- 1:ntr
  ttrans[ttrans==0] <- NA
  trans <- t(ttrans)
  return(trans)
}

