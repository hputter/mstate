crprep <- function(Tstop, status, data, trans=1, cens=0, Tstart=0, id, keep, shorten=TRUE, rm.na=TRUE, origin=0, prec.factor=100) {

  # Extract Tstop data
  if (!(is.numeric(Tstop))) {
    if (!is.character(Tstop) | length(Tstop)!=1)
      stop("argument \"Tstop\" should be a numeric vector or a character string")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"Tstop\" argument is a character string")
    tcol <- match(Tstop, names(data))
    if (is.na(tcol))
      stop("\"Tstop\" not found in data")
    Tstop <- data[, tcol]
  }
  else {
    if (!is.vector(Tstop))
        stop("argument should be a numeric vector or a character string")
  }

  nn <- length(Tstop)

  # Extract status data
  if (!(is.vector(status) & length(status) > 1)) {
    if (!is.character(status))
      stop("argument \"status\" should be a vector or a character string")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"status\" argument is a character string")
    tcol <- match(status, names(data))
    if (is.na(tcol))
      stop("\"status\" not found in data")
    status <- data[, tcol]
  }
  if (length(status) != nn)
    stop("Tstop and status have different lengths")

  # Extract Tstart data
  if (is.numeric(Tstart)&length(Tstart) == 1) Tstart <- rep(Tstart, nn)
  if (!(is.numeric(Tstart))) {
    if (!is.character(Tstart) | length(Tstart)!=1)
      stop("argument \"Tstart\" should be a numeric vector or a character string")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"Tstart\" argument is a character string")
    tcol <- match(Tstart, names(data))
    if (is.na(tcol))
      stop("\"Tstart\" not found in data")
    Tstart <- data[, tcol]
  }
  else {
    if (!is.vector(Tstart))
        stop("argument should be a numeric vector or a character string")
  }
  if (length(Tstart) != nn)
    stop("Tstop and Tstart have different lengths")

  sel <- !is.na(Tstart) & !is.na(Tstop)
  if (any(Tstart[sel] >= Tstop[sel]))
    stop("Tstop must be greater than Tstart")

  # Extract id data
  if (missing(id)) id <- patid <- 1:nn
  if (!(is.vector(id) & length(id) > 1)) {
    if (!is.character(id))
      stop("argument \"id\" should be a character vector")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"id\" argument is a character vector")
    tcol <- match(id, names(data))
    if (is.na(tcol))
      stop("\"id\" not found in data")
    patid <- data[, tcol]
    id <- 1:nn
  }
  else {
      if (length(id) != nn)
      stop("Tstop and id have different lengths")
      patid <- id
      id <- 1:nn
  }

  # Extract covariate data
  if(!missing(keep)) {
    if (!(is.matrix(keep) | is.data.frame(keep))) {
      if (is.character(keep)) {
        if (missing(data))
          stop("missing \"data\" argument not allowed when \"keep\" argument is a character vector")
        nkeep <- length(keep)
        kcols <- match(keep, names(data))
        if (any(is.na(kcols)))
          stop("at least one element of \"keep\" not found in data")
          keepname <- keep
          keep <- data[, kcols]
        }
        else {
          nkeep <- 1
          keepname <- names(keep)
          if (length(keep) != nn)
            stop("argument \"keep\" has incorrect dimensions")
        }
      }
      else {
        nkeep <- ncol(keep)
        if(is.data.frame(keep)){
            keepname <- names(keep)
        }
        else {
            keepname <- colnames(keep)
            if(is.null(keepname)) keepname <- paste("V",1:nkeep,sep="")
        }
        if (nrow(keep) != nn)
          stop("argument \"keep\" has incorrect dimensions")
        if (nkeep == 1)
          keep <- keep[, 1]
      }
  }

  Tstart <- Tstart - origin
  Tstop <- Tstop - origin

  # Eliminate records with missings in time or status
  if(rm.na){
    sel2 <- sel & !is.na(status)
    Tstart <- Tstart[sel2]
    Tstop <- Tstop[sel2]
    status <- status[sel2]
    id <- id[sel2]
    patid <- patid[sel2]
    n <- length(Tstop)
  }
  else {
    n <- nn
    sel2 <- sel
  }

  prec <- .Machine$double.eps*prec.factor
  calc.trunc <- any(Tstart != 0)
  # time.surv <- unique(sort(Tstop[status!=cens])) # unique event times of any type

  # Calculate time to entry (left truncation) distribution at t-
  if(calc.trunc) {
      surv.trunc <- survival::survfit(Surv(-Tstop,-Tstart+prec,rep(1,n))~1) # PL estimator of -L, corrected for ties wrt L and X (assume X > L)
      trunc.dist <- summary(surv.trunc)
      trunc.dist$time <- rev(-trunc.dist$time)
      trunc.dist$surv <- c(rev(trunc.dist$surv)[-1],1)
  }

  # Calculate product-limit time to censoring survival distribution, evaluated just before event time in case of ties
  surv.cens <- survival::survfit(Surv(Tstart,Tstop+ifelse(status==cens,prec,0),status==cens)~1) # PL censoring estimator
  cens.weights <- summary(surv.cens)

  # Create weighted data set for each event type of interest
  data.out <- list()
  i.list <- 1
  for(failcode in unique(trans)) {
    # unique event times of type of interest
    time.surv.T1 <- unique(sort(Tstop[status==failcode]))

    # Number of rows in new dataset
    Nw <- rep(1,n)
    tmp.compstat <- (status!=cens)&(status!=failcode)&!is.na(status)
    Nw[tmp.compstat] <-  apply(outer(Tstop[tmp.compstat],time.surv.T1,"<"),1,sum)+1

    # Create weighted dataset
    data.weight <- data.frame(id=rep(id,Nw), Tstart=NA, Tstop=NA, status=rep(status,Nw), weight.cens=NA)
    data.weight$Tstart <- unlist(lapply(1:n,FUN=function(x,tms,N) {if (N[x]==1) Tstart[x]
                                                                     else {if (N[x]==2) c(Tstart[x],Tstop[x]) else
                                                                     c(Tstart[x], Tstop[x], rev(rev(tms)[2:(N[x]-1)]))}}, tms=time.surv.T1,N=Nw))
    data.weight$Tstop <-  unlist(lapply(1:n,FUN=function(x,tms,N) {if (N[x]==1) Tstop[x]
                                                                     else c(Tstop[x], tail(tms,N[x]-1))}, tms=time.surv.T1, N=Nw))

    # add censoring weights to the new dataset
    data.weight$weight.cens <- unlist(lapply( split(data.weight$Tstop, data.weight$id),
                                      FUN=function(x, tms, srvs) if(length(x)==1) 1 else
                                        approx(tms, srvs, xout=x, method="constant", f=0, rule=2)$y/
                                        approx(tms, srvs, xout=x[1], method="constant", f=0, rule=2)$y, tms=c(0,cens.weights$time), srvs=c(1,cens.weights$surv)))

    # Add truncation weights to the new dataset
    if(calc.trunc) {
      weight.trunc <- unlist(lapply(split(data.weight$Tstop,data.weight$id),
                                         FUN=function(x,tms,srvs) if(length(x)==1) 1 else
                                           approx(tms, srvs, xout=x,method="constant",f=0,rule=2)$y/approx(tms, srvs,
                                           xout=x[1],method="constant",f=0,rule=2)$y,tms=trunc.dist$time,srvs=trunc.dist$surv))

      data.weight <- cbind(data.weight, weight.trunc)
    }

    tbl <- rle(data.weight$id)$length

    # Add covariates
    if(!missing(keep)) {

      # Extract covariate name from function
      if (is.null(keepname)) {
        m <- match.call(expand.dots = FALSE)
        m <- m[match("keep", names(m))]

        if(!is.null(m)) {
          keepname <- as.character(m[1])
          keepname.split <- strsplit(keepname, '')[[1]]

          tag <- which(keepname.split == '$')
          if(length(tag) != 0) {
            keepname <- substring(keepname, tag[length(tag)]+1)
          } else {
            tag <- which(keepname.split == '"')
            if(length(tag) != 0) {
              keepname <- substring(keepname, tag[1]+1, tag[2]-1)
            }
          }
        }
      }

      # Add covariates to the resultset
      if (nkeep > 0) {
        if (nkeep == 1) {
          keep <- keep[sel2]
          ddcovs <- rep(keep, tbl)
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- as.character(keepname)
        } else {
          keep <- keep[sel2, ]
          ddcovs <- lapply(1:nkeep, function(i) rep(keep[, i], tbl))
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- keepname
        }
        data.weight <- cbind(data.weight, ddcovs)
      }

    }

    # Shorten data set by combining rows with event types without censoring or truncation time in between
    if (shorten) {
      if(calc.trunc) {
        keep.rows <- c(diff(data.weight$id)!=0 | diff(data.weight$weight.cens)!=0 | diff(data.weight$weight.trunc)!=0, TRUE)
      } else {
        keep.rows <- c(diff(data.weight$id)!=0 | diff(data.weight$weight.cens)!=0, TRUE)
      }

      keep.rows[unlist(mapply(seq,1,tbl))==1] <- TRUE # First record always included in order to allow for CSH analysis
      keep.start <- data.weight$Tstart[unlist(lapply( split(keep.rows, data.weight$id), FUN=function(x) if(length(x)==1) x else c(TRUE,x[-length(x)])))]
      data.weight <- data.weight[keep.rows,]
      data.weight$Tstart <- keep.start
    }

    # Add count
    data.weight$count <- unlist(mapply(seq,1,rle(data.weight$id)$length)) # note: cannot use tbl if shorten=TRUE
    data.weight$failcode <- failcode
    if (shorten) {
         data.weight$id <- rep(patid,Nw)[keep.rows]
         } else
         data.weight$id <- rep(patid,Nw)

    data.out[[i.list]] <- data.weight
    i.list <- i.list+1
  }

  out <- do.call("rbind", data.out)
  row.names(out) <- as.character(1:nrow(out))
  class(out) <- c("crprep","data.frame")
  return(out)
}

