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

trans2tra <- function(trans)
  return(!(is.na(trans)))

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

tra2trans <- function(tra)
{
  ttrans <- t(matrix(as.numeric(tra), nrow(tra), ncol(tra)))
  ntr <- sum(tra)
  ttrans[ttrans==1] <- 1:ntr
  ttrans[ttrans==0] <- NA
  trans <- t(ttrans)
  return(trans)
}

