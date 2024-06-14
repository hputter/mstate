#' Function to prepare dataset for multi-state modeling in long format from
#' dataset in wide format
#' 
#' This function converts a dataset which is in wide format (one subject per
#' line, multiple columns indicating time and status for different states) into
#' a dataset in long format (one line for each transition for which a subject
#' is at risk). Selected covariates are replicated per subjects.
#' 
#' For \code{msprep}, the transition matrix should correspond to an
#' irreversible acyclic Markov chain. In particular, on the diagonals only
#' \code{NA}s are allowed.
#' 
#' The transition matrix, if irreversible and acyclic, will have starting
#' states, i.e. states into which no transitions are possible. For these
#' starting states \code{NA}s are allowed in the \code{time} and \code{status}
#' arguments, either as columns, when specified as matrix or data frame, or as
#' elements of the character vector when specified as character vector.
#' 
#' The function \code{msprep} uses a recursive algorithm through calls to the
#' recursive function \code{msprepEngine}. First, with the current transition
#' matrix, all starting states are detected (defined as states into which there
#' are no transitions). For each of these starting states, all subjects
#' starting from that state are selected and for each subject the next visited
#' state is detected by looking at all transitions from that starting state and
#' determining the smallest transition time with \code{status}=1. The recursive
#' \code{msprepEngine} is called again with the starting states deleted from
#' the transition matrix and with subjects deleted that either reached an
#' absorbing state or that were censored. For the remaining subjects the
#' starting states and times are updated in the next call. Datasets returned
#' from the \code{msprepEngine} calls are appended to the current dataset in
#' long format and finally sorted.
#' 
#' A warning is issued for a subject, if multiple transitions exist with the
#' same smallest transition time (and \code{status}=0). In such cases the next
#' transition cannot be determined unambiguously, and the state with the
#' smallest number is chosen. In our experience, occasionally the shortest
#' transition time has \code{status}=0, while a higher transition time has
#' \code{status}=1. Then this larger transition time and the corresponding
#' transition is selected. No warning is issued for these data inconsistencies.
#' 
#' @param time Either 1) a matrix or data frame of dimension n x S (n being the
#' number of individuals and S the number of states in the multi-state model),
#' containing the times at which the states are visited or last follow-up time,
#' or 2) a character vector of length S containing the column names indicating
#' these times. In the latter cases, some elements of \code{time} may be NA,
#' see Details
#' @param status Either 1) a matrix or data frame of dimension n x S,
#' containing, for each of the states, event indicators taking the value 1 if
#' the state is visited or 0 if it is not (censored), or 2) a character vector
#' of length S containing the column names indicating these status variables.
#' In the latter cases, some elements of \code{status} may be NA, see Details
#' @param data Data frame (not a tibble) in wide format in which to interpret
#' \code{time}, \code{status}, \code{id} or \code{keep}, if appropriate
#' @param trans Transition matrix describing the states and transitions in the
#' multi-state model. If S is the number of states in the multi-state model,
#' \code{trans} should be an S x S matrix, with (i,j)-element a positive
#' integer if a transition from i to j is possible in the multi-state model,
#' \code{NA} otherwise. In particular, all diagonal elements should be
#' \code{NA}. The integers indicating the possible transitions in the
#' multi-state model should be sequentially numbered, 1,...,K, with K the
#' number of transitions
#' @param start List with elements \code{state} and \code{time}, containing
#' starting states and times of the subjects in the data.  Default is
#' \code{NULL}, in which case all subjects start in state 1 at time 0. If a
#' single state and time are given this state and time is used for all
#' subjects, otherwise the length of \code{state} and \code{time} should equal
#' the number of subjects in \code{data}
#' @param id Either 1) a vector of length n containing the subject
#' identifications, or 2) a character string indicating the column name
#' containing these subject ids. If not provided, \code{"id"} will be assigned
#' with values 1,...,n
#' @param keep Either 1) a data frame or matrix with n rows or a numeric or
#' factor vector of length n containing covariate(s) that need to be retained
#' in the output dataset, or 2) a character vector containing the column names
#' of these covariates in \code{data}
#' @return An object of class \code{"msdata"}, which is a data frame in long
#' (counting process) format containing the subject id, the covariates
#' (replicated per subject), and \item{from}{the starting state} \item{to}{the
#' receiving state} \item{trans}{the transition number} \item{Tstart}{the
#' starting time of the transition} \item{Tstop}{the stopping time of the
#' transition} \item{status}{status variable, with 1 indicating an event
#' (transition), 0 a censoring} The \code{"msdata"} object has the transition
#' matrix as \code{"trans"} attribute.
#' @author Hein Putter \email{H.Putter@@lumc.nl} and Marta Fiocco
#' @references Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics:
#' Competing risks and multi-state models. \emph{Statistics in Medicine}
#' \bold{26}, 2389--2430.
#' @keywords datagen
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # some data in wide format
#' tg <- data.frame(stt=rep(0,6),sts=rep(0,6),
#'         illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'         dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'         x1=c(1,1,1,2,2,2),x2=c(6:1))
#' tg$x1 <- factor(tg$x1,labels=c("male","female"))
#' tg$patid <- factor(2:7,levels=1:8,labels=as.character(1:8))
#' # define time, status and covariates also as matrices
#' tt <- matrix(c(rep(NA,6),tg$illt,tg$dt),6,3)
#' st <- matrix(c(rep(NA,6),tg$ills,tg$ds),6,3)
#' keepmat <- data.frame(gender=tg$x1,age=tg$x2)
#' # data in long format using msprep
#' msprep(time=tt,status=st,trans=tmat,keep=as.matrix(keepmat))
#' msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),data=tg,
#' 		id="patid",keep=c("x1","x2"),trans=tmat)
#' # Patient no 5, 6 now start in state 2 at time t=4 and t=10
#' msprep(time=tt,status=st,trans=tmat,keep=keepmat,
#'         start=list(state=c(1,1,1,1,2,2),time=c(0,0,0,0,4,10)))
#' 
#' @export msprep
msprep <- function (time, status, data, trans, start, id, keep) 
{
    if (is.circular(trans))
        stop("msprep cannot be used with a circular transition matrix")
    if (!missing(data)) data <- as.data.frame(data)
    if (!(is.matrix(time) | (is.data.frame(time)))) {
        if (!is.character(time)) 
            stop("argument \"time\" should be a character vector")
        if (missing(data)) 
            stop("missing \"data\" argument not allowed when time argument is character vector")
        startings <- which(apply(!is.na(trans), 2, sum) == 0)
        wh <- which(is.na(time))
        if (!all(wh %in% startings)) 
            stop("no NA's allowed in the \"time\" argument for non-starting states")
        tcols <- match(time[!is.na(time)], names(data))
        if (any(is.na(tcols))) 
            stop("at least one of elements of \"time\" not in data")
        time <- matrix(NA, nrow(data), length(time))
        whcols <- (1:ncol(time))[!(1:ncol(time)) %in% wh]
        time[, whcols] <- as.matrix(data[, tcols])
    }
    if (!(is.matrix(status) | (is.data.frame(status)))) {
        if (!is.character(status)) 
            stop("argument \"status\" should be a character vector")
        if (missing(data)) 
            stop("missing \"data\" argument not allowed when status argument is character vector")
        startings <- which(apply(!is.na(trans), 2, sum) == 0)
        wh <- which(is.na(status))
        if (!all(wh %in% startings)) 
            stop("no NA's allowed in the \"status\" argument for non-starting states")
        dcols <- match(status[!is.na(status)], names(data))
        if (any(is.na(dcols))) 
            stop("at least one of elements of \"status\" not in data")
        status <- matrix(NA, nrow(data), length(status))
        whcols <- (1:ncol(status))[!(1:ncol(status)) %in% wh]
        status[, whcols] <- as.matrix(data[, dcols])
    }
    time <- as.matrix(time)
    status <- as.matrix(status)
    if (!all(dim(time) == dim(status))) 
        stop("unequal dimensions of \"time\" and \"status\" data")
    n <- nrow(time)
    K <- dim(trans)[1]
    if ((dim(trans)[2] != K) | (dim(time)[2] != K)) 
        stop("dimension of \"trans\" does not match with length of \"time\" and \"status\"")
    idname <- "id"
    missingid <- FALSE
    if (missing(id)) {
        missingid <- TRUE
        id <- 1:n
    }
    else {
        if (!is.vector(id)) 
            stop("argument \"id\" is not a vector")
        else {
            if (!is.character(id)) {
                if (length(id) != n) 
                  stop("argument \"id\" not of correct length")
            }
            else {
                if (length(id) == 1) {
                  if (n == 1) 
                    stop("cannot determine whether \"id\" argument indicates ")
                  else {
                    idname <- id
                    id <- as.factor(data[[id]])
                  }
                }
                else {
                  if (length(id) != n) 
                    stop("argument \"id\" not of correct length")
                  id <- factor(id)
                }
            }
        }
    }
    idlevels <- NULL
    if (is.factor(id)) 
        idlevels <- levels(id)
    if (!missing(start)) {
        startstate <- start$state
        starttime <- start$time
        if (length(startstate) != length(starttime)) 
            stop("starting states and times not of equal length")
        if (length(startstate) > 1) {
            if (length(startstate) != n) 
                stop("length of starting states and times different from no of subjects in data")
        }
        else {
            startstate <- rep(startstate, n)
            starttime <- rep(starttime, n)
        }
    }
    else {
        startstate <- rep(1, n)
        starttime <- rep(0, n)
    }
    
    ord1 <- order(id)
    
    msres <- msprepEngine(time = time, status = status, id = id, 
        starttime = starttime, startstate = startstate, trans = trans, 
        originalStates = (1:nrow(trans)), longmat = NULL)
    msres <- as.data.frame(msres)
    names(msres) <- c(idname, "from", "to", "trans", "Tstart", 
        "Tstop", "status")
    msres$time <- msres$Tstop - msres$Tstart
    msres <- msres[, c(1:6, 8, 7)]
    ord <- order(msres[, 1], msres[, 5], msres[, 2], msres[, 3])
    msres <- msres[ord, ]
    row.names(msres) <- 1:nrow(msres)
    if (!is.null(idlevels)) 
        msres[, 1] <- factor(as.integer(msres[, 1]), 1:length(idlevels), 
            labels = idlevels)
    if (!missing(keep)) {
        if (!(is.matrix(keep) | (is.data.frame(keep)))) {
            if (is.character(keep)) {
                if (missing(data)) 
                  stop("argument \"data\" is missing, with no default")
                nkeep <- length(keep)
                kcols <- match(keep, names(data))
                if (any(is.na(kcols))) 
                  stop("at least one of elements of \"keep\" not in data")
                keepname <- keep
                keep <- data[, kcols]
            }
            else {
                nkeep <- 1
                keepname <- names(keep)
                if (is.null(keepname)) keepname <- "keep"
                if (length(keep) != n) 
                  stop("argument \"keep\" has incorrect dimension")
            }
        }
        else {
            nkeep <- ncol(keep)
            keepname <- names(keep)
            if (nrow(keep) != n) 
                stop("argument \"keep\" has incorrect dimension")
            if (nkeep == 1) 
                keep <- keep[, 1]
        }
        if (is.null(keepname))
            keepname <- paste("keep", as.character(1:nkeep), 
                sep = "")
        if (nkeep > 0) {
            if (is.factor(msres[, 1])) 
                msres[, 1] <- factor(msres[, 1])
            tbl <- table(msres[, 1])
            if (nkeep > 1)
              keep <- keep[ord1,,drop=FALSE]
            if (nkeep == 1) {
              ddcovs <- rep(keep[ord1], tbl)
              ddcovs <- as.data.frame(ddcovs)
              names(ddcovs) <- keepname
            }
            else {
              ddcovs <- lapply(1:nkeep, function(i) rep(keep[, 
                i], tbl))
              ddcovs <- as.data.frame(ddcovs)
              names(ddcovs) <- keepname
            }
            msres <- cbind(msres, ddcovs)
        }
    }
    attr(msres, "trans") <- trans
    class(msres) <- c("msdata", "data.frame")
    return(msres)
}

msprepEngine <- function(time,status,id,starttime,startstate,trans,originalStates,longmat)
{
### Recursive engine for msprep
### Input:
###       time: n x K numeric matrx containing arrival times
###       status: n x K numeric matrix containing arrival times
###       id: n vector with ids
###       starttime: n vector of starting times
###       startstate: n vector of starting states
###       trans: current (K x K) matrix
###       originalStates: the numbers of the original states in the
###         current transition matrix
###       longmat: the matrix in longformat that has already been
###         constructed
### Output:
###       A new longmat matrix with new data appended 
    if (is.null(nrow(time))) return(longmat) # finished
    if (nrow(time)==0) return(longmat) # also finished
    states.to <- apply(!is.na(trans),1,sum)
    absorbing <- which(states.to==0)
    states.from <- apply(!is.na(trans),2,sum)
    startings <- which(states.from==0)
    # preset values for next call to msprepEngine
    newstate <- startstate
    newtime <- starttime
    to.remove <- NULL # indices in data, state, time to be removed at the end
    for (starting in startings) {
        # select all subjects starting in starting, no is nstart
        subjs <- which(startstate==starting)
        nstart <- length(subjs)
        # determine states that can be reached from starting
        tostates <- which(!is.na(trans[starting,]))
        transs <- trans[starting,tostates]
        nreach <- length(tostates)
        # make matrix with nstart*nreach rows
        if ((nstart>0) & (nreach>0)) {
            Tstart <- starttime[subjs]
            Tstop <- time[subjs,tostates,drop=FALSE]
            # event times before start time do not count
            # the statement below defining hlp makes sure
            # these are considered as censorings
            Tstop[Tstop<Tstart] <- Inf
            stat <- status[subjs,tostates,drop=FALSE]
            smallesttime <- apply(Tstop,1,min) # determine first event time
            # now determine which is the corresponding state
            # bit tricky because of censoring and possible data errors
            # limited warnings given at present ...
            hlp <- Tstop * 1/stat
            hlp[Tstop==0 & stat==0] <- Inf # define 0/0 as infinity
            nexttime <- apply(hlp,1,min) # determine first event time
            censored <- which(is.infinite(apply(hlp,1,min)))
            wh <- which(smallesttime<nexttime)
            whminc <- setdiff(wh,censored)
            if (length(whminc)>0) {
                whsubjs <- id[subjs[whminc]]; whsubjs <- paste(whsubjs,collapse=" ")
                warning("From starting state ",originalStates[starting],", subject ",
                    whsubjs," has smallest transition time with status=0, larger transition time with status=1")
            }
            nexttime[censored] <- smallesttime[censored]
            if (ncol(hlp)>1) {
                hlpsrt <- t(apply(hlp,1,sort))
                warn1 <- which(hlpsrt[,1]-hlpsrt[,2]==0)
                if (length(warn1)>0)
                {
                    isw <- id[subjs[warn1]]; isw <- paste(isw,collapse=" ")
                    hsw <- hlpsrt[warn1,1]; hsw <- paste(hsw,collapse=" ")
                    warning("Starting from state ",originalStates[starting],
                        ", simultaneous transitions possible for subjects ",
                        isw," at times ",hsw,
                        "; smallest receiving state chosen")
                }
            }
            if (length(censored)>0) {
                nextstate <- apply(hlp[-censored,,drop=FALSE],1,which.min)
                reachAbsorb <- (1:nstart)[-censored][which(tostates[nextstate] %in% absorbing)]
            } else {
                nextstate <- apply(hlp,1,which.min)
                reachAbsorb <- (1:nstart)[which(tostates[nextstate] %in% absorbing)]
            }
            # the status to be returned in long dataframe has 0 if censored
            # and 1 for the transition followed
            statmat <- matrix(0,nstart,nreach)
            if (length(censored)>0) statmatmin <- statmat[-censored,,drop=FALSE] else statmatmin <- statmat
            if (nrow(statmatmin)>0)
                statmatmin <- t(sapply(1:nrow(statmatmin),function(i) {
                    x <- statmatmin[i,]
                    x[nextstate[i]] <- 1
                    return(x) }
                ))
            if (length(censored)>0) statmat[-censored,] <- statmatmin else statmat <- statmatmin
            mm <- matrix(c(
                rep(id[subjs],rep(nreach,nstart)), # id
                rep(originalStates[starting],nreach*nstart), # from
                rep(originalStates[tostates],nstart), # to
                rep(transs,nstart), # trans
                rep(Tstart,rep(nreach,nstart)), # Tstart
                rep(nexttime,rep(nreach,nstart)), # Tstop
                as.vector(t(statmat))), # status
                nreach*nstart,7)
            # stack upon what is already there
            longmat <- rbind(longmat,mm)
            # adjust data: remove subjects who didn't reach any new state
            #              for those who do reach new state, adjust starting
            #              state and time
            to.remove <- c(to.remove,subjs[c(censored,reachAbsorb)])
            if (length(censored)>0) newstate[subjs[-censored]] <- tostates[nextstate]
            else newstate[subjs] <- tostates[nextstate]
            if (length(censored)>0) newtime[subjs[-censored]] <- nexttime[-censored]
            else newtime[subjs] <- nexttime
        }
    }
    if (length(to.remove)>0) {
        time <- time[-to.remove, , drop = FALSE]
        status <- status[-to.remove, , drop = FALSE]
        newtime <- newtime[-to.remove]          
        newstate <- newstate[-to.remove]
        id <- id[-to.remove]
    }
    # Some states will be removed from the transition matrix in the next call
    # We have to adjust newstate to mean the new states
    K <- nrow(trans)
    idx <- rep(1,K); idx[startings] <- 0; idx <- cumsum(idx)
    newstate <- idx[newstate]
    Recall(time=time[,-startings],status=status[,-startings],
            id=id,starttime=newtime,startstate=newstate,
            trans=trans[-startings,-startings],
            originalStates=originalStates[-startings],
            longmat=longmat)
}
