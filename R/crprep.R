#' Function to create weighted data set for competing risks analyses
#' 
#' This function converts a dataset that is in short format (one subject per
#' line) into a counting process format with time-varying weights that correct
#' for right censored and left truncated data. With this data set, analyses
#' based on the subdistribution hazard can be performed.
#' 
#' For each event type as specified via \code{trans}, individuals with a
#' competing event remain in the risk set with weights that are determined by
#' the product-limit forms of the time-to-censoring and time-to-entry
#' estimates. Typically, their weights change over follow-up, and therefore
#' such individuals are split into several rows. Censoring weights are always
#' computed. Truncation weights are computed only if \code{Tstart} is
#' specified.
#' 
#' If several event types are specified at once, regression analyses using the
#' stacked format data set can be performed (see Putter et al. 2007 and Chapter
#' 4 in Geskus 2016). The data set can also be used for a regression on the
#' cause-specific hazard by restricting to the subset \code{subset=count==0}.
#' 
#' Missing values are allowed in \code{Tstop}, \code{status}, \code{Tstart},
#' \code{strata} and \code{keep}. Rows for which \code{Tstart} or \code{Tstart}
#' is missing are deleted.
#' 
#' There are two ways to supply the data. If given "by value" (option 1), the
#' actual data vectors are used. If given "by name" (option 2), the column
#' names are specified, which are read from the data set in \code{data}. In
#' general, the second option is preferred.
#' 
#' If data are given by value, the following holds for the naming of the
#' columns in the output data set. If \code{keep}, \code{strata} or \code{id}
#' is a vector from a (sub)-list, e.g. obj$name2$name1, then the column name is
#' based on the most inner part (i.e.\ "name1"). If it is a vector of the form
#' obj[,"name1"], then the column is named "name1". For all other vector
#' specifications, the name is copied as is. If \code{keep} is a data.frame or
#' a named matrix, the same names are used for the covariate columns in the
#' output data set. If keep is a matrix without names, then the covariate
#' columns are given the names "V1" until "Vk".
#' 
#' The current function does not allow to create a weighted data set in which
#' the censoring and/or truncation mechanisms depend on covariates via a
#' regression model.
#' 
#' @aliases crprep crprep.default
#' @param Tstop Either 1) a vector containing the time at which the follow-up
#' is ended, or 2) a character string indicating the column name in \code{data}
#' that contains the end times (see Details).
#' @param status Either 1) a vector describing status at end of follow-up,
#' having the same length as \code{Tstop}, or 2) a character string indicating
#' the column name that contains this information.
#' @param data Data frame in which to interpret \code{Tstart}, \code{status},
#' \code{Tstart}, \code{id}, \code{strata} and \code{keep}, if given as
#' character value (specification 2, "by name").
#' @param trans Values of \code{status} for which weights are to be calculated.
#' @param cens Value that denotes censoring in \code{status} column.
#' @param Tstart Either 1) a vector containing the time at which the follow-up
#' is started, having the same length as \code{Tstop}, or 2) a character string
#' indicating the column name that contains the entry times, or 3) one numeric
#' value in case it is the same for every subject. Default is 0.
#' @param id Either 1) a vector, having the same length as \code{Tstop},
#' containing the subject identifiers, or 2) a character string indicating the
#' column name containing these subject identifiers. If not provided, a column
#' \code{id} is created with subjects having values 1,...,n.
#' @param strata Either 1) a vector of the same length as \code{Tstop}, or 2) a
#' character string indicating the column name that contains this information.
#' Weights are calculated for per value in this vector.
#' @param keep Either 1) a data frame or matrix or a numeric or factor vector
#' containing covariate(s) that need to be retained in the output dataset.
#' Number of rows/length should correspond with \code{Tstop}, or 2) a character
#' vector containing the column names of these covariates in \code{data}.
#' @param shorten Logical. If true, number of rows in output is reduced by
#' collapsing rows within a subject in which weights do not change.
#' @param rm.na Logical. If true, rows for which \code{status} is missing are
#' deleted.
#' @param origin Substract origin time units from all Tstop and Tstart times.
#' @param prec.factor Factor by which to multiply the machine's precision.
#' Censoring and truncation times are shifted by prec.factor*precision if event
#' times and censoring/truncation times are equal.
#' @param list() further arguments to be passed to or from other methods. They
#' are ignored in this function.
#' @return A data frame in long (counting process) format containing the
#' covariates (replicated per subject). The following column names are used:
#' \item{Tstart}{start dates of dataset} \item{Tstop}{stop dates of dataset}
#' \item{status}{status of the subject at the end of that row}
#' \item{weight.cens}{weights due to censoring mechanism}
#' \item{weight.trunc}{weights due to truncation mechanism (if present)}
#' \item{count}{row number within subject and event type under consideration}
#' \item{failcode}{event type under consideration}
#' 
#' The first column is the subject identifier. If the argument "id" is missing,
#' it has values 1:n and is named "id". Otherwise the information is taken from
#' the \code{id} argument.
#' 
#' Variables as specified in \code{strata} and/or \code{keep} are included as
#' well (see Details).
#' @author Ronald Geskus
#' @references Geskus RB (2011). Cause-Specific Cumulative Incidence Estimation
#' and the Fine and Gray Model Under Both Left Truncation and Right Censoring.
#' \emph{Biometrics} \bold{67}, 39--49.
#' 
#' Geskus, Ronald B. (2016). \emph{Data Analysis with Competing Risks and
#' Intermediate States.} CRC Press, Boca Raton.
#' 
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: Competing
#' risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
#' 2389--2430.
#' @keywords datagen survival
#' @examples
#' 
#' data(aidssi)
#' aidssi.w <- crprep("time", "cause", data=aidssi, trans=c("AIDS","SI"),
#'                    cens="event-free", id="patnr", keep="ccr5")
#' 
#' # calculate cause-specific cumulative incidence, no truncation,
#' # compare with Cuminc (also from mstate)
#' ci <- Cuminc(aidssi$time, aidssi$status)
#' sf <- survfit(Surv(Tstart,Tstop,status=="AIDS")~1, data=aidssi.w,
#'               weight=weight.cens, subset=failcode=="AIDS")
#' plot(sf, fun="event", mark.time=FALSE)
#' lines(CI.1~time,data=ci,type="s",col="red")
#' sf <- survfit(Surv(Tstart,Tstop,status=="SI")~1, data=aidssi.w,
#'               weight=weight.cens, subset=failcode=="SI")
#' plot(sf, fun="event", mark.time=FALSE)
#' lines(CI.2~time,data=ci,type="s",col="red")
#' 
#' # Fine and Gray regression for cause 1
#' cw <- coxph(Surv(Tstart,Tstop,status=="AIDS")~ccr5, data=aidssi.w,
#'       weight=weight.cens, subset=failcode=="AIDS")
#' cw
#' # This can be checked with the results of crr (cmprsk)
#' # crr(ftime=aidssi$time, fstatus=aidssi$status, cov1=as.numeric(aidssi$ccr5))
#' 
#' # Gray's log-rank test
#' aidssi.wCCR <- crprep("time", "cause", data=aidssi, trans=c("AIDS","SI"),
#'                       cens="event-free", id="patnr", strata="ccr5")
#' test.AIDS <- coxph(Surv(Tstart,Tstop,status=="AIDS")~ccr5, data=aidssi.wCCR,
#'                    weights=weight.cens, subset=failcode=="AIDS")
#' test.SI <- coxph(Surv(Tstart,Tstop,status=="SI")~ccr5, data=aidssi.wCCR,
#'                  weights=weight.cens, subset=failcode=="SI")
#' ## score test statistic and p-value
#' c(test.AIDS$score, 1-pchisq(test.AIDS$score,1)) # AIDS
#' c(test.SI$score, 1-pchisq(test.SI$score,1))     # SI
#' # This can be compared with the results of cuminc (cmprsk)
#' # with(aidssi, cuminc(time, status, group=ccr5)$Tests)
#' # Note: results are not exactly the same
#' 
#' @export 
crprep <- function(Tstop, ...) UseMethod("crprep")

#' @method crprep default
#' @export
crprep.default <-
function(Tstop, status, data, trans=1, cens=0, Tstart=0, id, strata, keep, shorten=TRUE, rm.na=TRUE, origin=0,
         prec.factor=1000, ...) {

  ## Extract Tstop data if given by column name
  if (!(is.numeric(Tstop))) {
    if (!is.character(Tstop) | length(Tstop)!=1)
      stop("argument \"Tstop\" should be a numeric vector or a character string")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"Tstop\" argument is a character string")
    tcol <- match(Tstop, names(data))
    if (is.na(tcol))
      stop("\"Tstop\" not found in data")
    Tstop <- data[, tcol]
  } else {
    if (!is.vector(Tstop))
      stop("argument should be a numeric vector or a character string")
  }
  nn <- length(Tstop)

  ## Extract Tstart data if given by column name
  if (is.numeric(Tstart)&length(Tstart) == 1) {
    Tstart <- rep(Tstart, nn)
  } else {
    if (!(is.numeric(Tstart))) {
      if (!is.character(Tstart) | length(Tstart)!=1)
        stop("argument \"Tstart\" should be a numeric vector or a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"Tstart\" argument is a character string")
      tcol <- match(Tstart, names(data))
      if (is.na(tcol))
        stop("\"Tstart\" not found in data")
      Tstart <- data[, tcol]
      } else {
      if (!is.vector(Tstart))
        stop("argument should be a numeric vector or a character string")
      }
  }
  if (length(Tstart) != nn)
    stop("Tstop and Tstart have different lengths")
  ## Check whether Tstart is needed
  calc.trunc <- any(Tstart[!is.na(Tstart)] != 0)

  ## Select rows without missing time value
  sel <- !is.na(Tstart) & !is.na(Tstop)
  if (any(Tstart[sel] >= Tstop[sel]))
    stop("Tstop must be greater than Tstart")

  ## Extract status data if given by column name
  if (length(status) == 1) {
    if (!is.character(status))
      stop("argument \"status\" should be a vector or a character string")
    if (missing(data))
      stop("missing \"data\" argument not allowed when \"status\" argument is a character string")
    tcol <- match(status, names(data))
    if (is.na(tcol))
      stop("\"status\" not found in data")
    status <- data[ ,tcol]
  }
  if (length(status) != nn)
    stop("Tstop and status have different lengths")

  ## Extract strata data; value 1 if not specified
  if (missing(strata)) {
    strata.val <- rep(1,nn)
  } else {
    if (is.matrix(strata) | is.data.frame(strata))
        stop("only one variable is allowed in \"strata\"")
    if (!(is.vector(as.numeric(factor(strata))) & length(strata) > 1)) {
      if (!is.character(strata))
        stop("argument \"strata\" should be a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"strata\" argument is a character string")
      tcol <- match(strata, names(data))
      if (is.na(tcol))
        stop("\"strata\" not found in data")
      strata.name <- strata
      strata.val <- data[ ,tcol]
    } else {
      if (length(strata) != nn)
        stop("Tstop and strata have different lengths")
      strata.name <- names(strata)
      strata.val <- strata
    }
  }
  strata.num <- as.numeric(factor(strata.val))

  ## Extract id data; values 1:nn if not specified
  if (missing(id)) {
    id.name <- "id"
    id  <- num.id <- 1:nn
  } else {
    if (is.matrix(id) | is.data.frame(id))
      stop("only one variable is allowed in \"id\"")
    if (!(is.vector(id) & length(id) > 1)) { # by name
      if (!is.character(id))
        stop("argument \"id\" should be a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"id\" argument is a character string")
      tcol <- match(id, names(data))
      if (is.na(tcol))
        stop("\"id\" not found in data")
      id.name <- id
      num.id <- 1:nn
      id <- data[, tcol]
    } else {                                 # by value
      if (length(id) != nn)
        stop("Tstop and id have different lengths")
      id.name <- names(id)
      num.id <- 1:nn
    }
  }

  ## Eliminate records with missings in status if rm.na=TRUE
  if(rm.na) sel <- sel & !is.na(status)
  Tstart <- Tstart[sel]
  Tstop <- Tstop[sel]
  status <- status[sel]
  strata.val <- strata.val[sel]
  strata.num <- strata.num[sel]
  id <- id[sel]
  num.id <- num.id[sel]
  n <- length(Tstop)

  ## Extract covariate data
  if(!missing(keep)) {
    if (!(is.matrix(keep) | is.data.frame(keep))) {
      if (is.character(keep)) {  # if given by column name
        if (missing(data))
          stop("missing \"data\" argument not allowed when \"keep\" argument is a character vector")
        nkeep <- length(keep)
        kcols <- match(keep, names(data))
        if (any(is.na(kcols)))
          stop("at least one element of \"keep\" not found in data")
        keep.name <- keep
        keep <- data[, kcols]
      } else {                     # if one column, given by value
        nkeep <- 1
##        keep.name <- names(as.data.frame(keep))
        keep.name <- names(keep)
##        if(is.null(keep.name)) keep.name <- "V1"
        if (length(keep) != nn)
          stop("Tstop and keep have different lengths")
      }
    } else {                         # if a matrix/data.frame
      nkeep <- ncol(keep)
      if(is.data.frame(keep)){
        keep.name <- names(keep)
      } else {
        keep.name <- colnames(keep)
        if(is.null(keep.name)) keep.name <- paste("V",1:nkeep,sep="")
      }
      if (nrow(keep) != nn)
        stop("length Tstop and number of rows in keep are differents")
      if (nkeep == 1)
        keep <- keep[, 1]
    }
  }

  Tstart <- Tstart - origin
  Tstop <- Tstop - origin


  ## Start calculations
  prec <- .Machine$double.eps*prec.factor

  ## Calculate product-limit time-to-censoring distribution, "event" not included in case of ties
  surv.cens <- survival::survfit(Surv(Tstart,Tstop+ifelse(status==cens,prec,0),status==cens)~strata.num)

  ## Calculate time to entry (left truncation) distribution at t-, use 2*prec in order to exclude censorings at same time
  if(calc.trunc) surv.trunc <- survival::survfit(Surv(-Tstop,-(Tstart+2*prec),rep(1,n))~strata.num)
  ## trunc.dist <- summary(surv.trunc)
  ## trunc.dist$time <- rev(-trunc.dist$time)-prec
  ## trunc.dist$surv <- c(rev(trunc.dist$surv)[-1],1)

  ## Create weighted data set for each event type as specified in trans
  data.out <- vector("list",length(trans))
  i.list <- 1
  strat <- sort(unique(strata.num),na.last=TRUE)
  len.strat <- length(strat)
  ## Start weight calculation per event type
  for(failcode in trans) {
    if(len.strat==1){ # no strata
      data.weight <- create.wData.omega(Tstart, Tstop, status, num.id, 1, failcode, cens)
      tmp.time <- data.weight$Tstop
      data.weight$weight.cens[order(tmp.time)] <- summary(surv.cens, times=tmp.time-prec)$surv
      if(calc.trunc) data.weight$weight.trunc[order(-tmp.time)] <- summary(surv.trunc, times=-tmp.time)$surv
    } else {
      data.weight <- vector("list",len.strat)
      if(is.na(strat[len.strat])) {
        tmp.sel <- is.na(strata.num)
        data.weight[[len.strat]] <- data.frame(id=num.id[tmp.sel], Tstart=Tstart[tmp.sel], Tstop=Tstop[tmp.sel], status=status[tmp.sel], strata=NA,  weight.cens=NA)
        if(calc.trunc) data.weight[[len.strat]]$weight.trunc <- NA
      }
      for(tmp.strat in 1:(len.strat-is.na(strat[len.strat]))){
        tmp.sel <- !is.na(strata.num) & strata.num==tmp.strat
        data.weight[[tmp.strat]] <- create.wData.omega(Tstart[tmp.sel], Tstop[tmp.sel], status[tmp.sel], num.id[tmp.sel], tmp.strat, failcode, cens)
        tmp.time <- data.weight[[tmp.strat]]$Tstop
        data.weight[[tmp.strat]]$weight.cens[order(tmp.time)] <- summary(surv.cens[tmp.strat], times=tmp.time-prec)$surv
        if(calc.trunc) data.weight[[tmp.strat]]$weight.trunc[order(-tmp.time)] <- summary(surv.trunc[tmp.strat], times=-tmp.time)$surv
      }
      data.weight <- do.call("rbind", data.weight)
    }
    ## Calculate omega-censoring weights
    data.weight <- data.weight[order(data.weight$id,data.weight$Tstop), ]
    data.weight$weight.cens <- unlist(tapply(data.weight$weight.cens, data.weight$id, FUN=function(x) if(length(x)==1&!is.na(x[1])) 1 else x/x[1]))

    ## Calculate omega-truncation weights
    if(calc.trunc) {
      data.weight$weight.trunc <- unlist(tapply(data.weight$weight.trunc, data.weight$id, FUN=function(x) if(length(x)==1&!is.na(x[1])) 1 else x/x[1]))
    }

    tbl <- table(data.weight$id)

    ## Add covariates
    if(!missing(keep)) {
      ## Extract covariate name from function
      if (is.null(keep.name)) {
        m <- match.call(expand.dots = FALSE)
        m <- m[match("keep", names(m))]
        if(!is.null(m)) {
          keep.name <- as.character(m[1])
          keep.name.split <- strsplit(keep.name, '')[[1]]
          tag <- which(keep.name.split == '$')
          if(length(tag) != 0) {
            keep.name <- substring(keep.name, tag[length(tag)]+1)
          } else {
            tag <- which(keep.name.split == '"')
            if(length(tag) != 0) {
              keep.name <- substring(keep.name, tag[1]+1, tag[2]-1)
            }
          }
        }
      }

      ## Add covariates to the resultset
      if (nkeep > 0) {
        if (nkeep == 1) {
          keep <- keep[sel]
          ddcovs <- rep(keep, tbl)
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- as.character(keep.name)
        } else {
          keep <- keep[sel, ]
          ddcovs <- lapply(1:nkeep, function(i) rep(keep[, i], tbl))
          ddcovs <- as.data.frame(ddcovs)
          names(ddcovs) <- keep.name
        }
        data.weight <- cbind(data.weight, ddcovs)
      }
    }

    ## Shorten data set by combining rows with event types without censoring or truncation time in between
    if (shorten) {
      if(calc.trunc) {
        keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0 | diff(weight.trunc)!=0, TRUE))
      } else {
        keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0, TRUE))
      }
      ## First record always included as separate row in order to allow for CSH analysis
      keep.rows[unlist(mapply(seq,1,tbl))==1] <- TRUE
      keep.start <- data.weight$Tstart[unlist(tapply(keep.rows, data.weight$id, FUN=function(x) if(length(x)==1) x else c(TRUE,x[-length(x)])))]
      data.weight <- data.weight[keep.rows,]
      data.weight$Tstart <- keep.start
    }

    ## Recalculate tbl after shorten
    tbl <- table(data.weight$id)
    ## Add count
    data.weight$count <- unlist(mapply(seq,1,tbl))
    data.weight$failcode <- failcode
    ## Return to original id
    data.weight$id <- rep(id,tbl)

    data.out[[i.list]] <- data.weight
    i.list <- i.list+1
  }

  out <- do.call("rbind", data.out)

  if(!missing(strata)) {
  ## Extract strata name if given by value
    if (is.null(strata.name)) {
      m <- match.call(expand.dots = FALSE)
      m <- m[match("strata", names(m))]
      if(!is.null(m)) {
        strata.name <- as.character(m[1])
        strata.name.split <- strsplit(strata.name, '')[[1]]
        tag <- which(strata.name.split == '$')
        if(length(tag) != 0) {
          strata.name <- substring(strata.name, tag[length(tag)]+1)
        } else {
          tag <- which(strata.name.split == '"')
          if(length(tag) != 0) {
            strata.name <- substring(strata.name, tag[1]+1, tag[2]-1)
          }
        }
      }
    }
    ## Use original stratum values
    if(is.factor(strata.val)) {
      out$strata <- factor(out$strata, labels=levels(strata.val))
    } else {
      out$strata <- levels(factor(strata.val))[out$strata]
      if(is.numeric(strata.val)) out$strata <- as.numeric(out$strata)
    }
    ## Use original name of column
    tmp.sel <- match("strata", names(out))
    names(out)[tmp.sel] <- strata.name
  } else {
    out$strata <- NULL
  }

  if (is.null(id.name)) {
    m <- match.call(expand.dots = FALSE)
    m <- m[match("id", names(m))]
    if(!is.null(m)) {
      id.name <- as.character(m[1])
      id.name.split <- strsplit(id.name, '')[[1]]
      tag <- which(id.name.split == '$')
      if(length(tag) != 0) {
        id.name <- substring(id.name, tag[length(tag)]+1)
      } else {
        tag <- which(id.name.split == '"')
        if(length(tag) != 0) {
          id.name <- substring(id.name, tag[1]+1, tag[2]-1)
        }
      }
    }
  }

  row.names(out) <- as.character(1:nrow(out))
  names(out)[1] <- id.name
  class(out) <- c("crprep","data.frame")
  return(out)
}

