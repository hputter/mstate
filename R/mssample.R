#' Sample paths through a multi-state model
#' 
#' Given cumulative transition hazards sample paths through the multi-state
#' model.
#' 
#' The procedure is described in detail in Fiocco, Putter & van Houwelingen
#' (2008). The argument \code{beta.state} and the element \code{tstate} from
#' the argument \code{history} are meant to incorporate situations where the
#' time at which some previous states were visited may affect future transition
#' rates. The relation between time of visit of state \code{s} and transition
#' \code{k} is assumed to be linear on the log-hazards; the corresponding
#' regression coefficient is to be supplied as the (s,k)-element of
#' \code{beta.state}, which is 0 if no such effect has been included in the
#' model. If no such effects are present, then \code{beta.state}=\code{NULL}
#' (default) suffices. In the \code{tstate} element of \code{history}, the
#' \code{s}-th element is the time at which state \code{s} was visited. This is
#' only relevant for states which have been visited prior to the beginning of
#' sampling, i.e. before the \code{time} element of \code{history}; the
#' elements of \code{tstate} are internally updated when in the sampling
#' process new states are visited (only if \code{beta.state} is not \code{NULL}
#' to avoid unnecessary computations).
#' 
#' @param Haz Cumulative hazards to be sampled from. These should be given as a
#' data frame with columns \code{time}, \code{Haz}, \code{trans}, for instance
#' as the \code{Haz} list element given by \code{\link{msfit}}.
#' @param trans Transition matrix describing the multi-state model.  See
#' \code{trans} in \code{\link{msprep}} for more detailed information
#' @param history A list with elements \code{state}, specifying the starting
#' state(s), \code{time}, the starting time(s), and \code{tstate}, a numeric
#' vector of length the number of states, specifying at what times states have
#' been visited, if appropriate. The default of \code{tstate} is \code{NULL};
#' more information can be found under Details.
#' 
#' The elements \code{state} and \code{time} may either be scalars or vectors,
#' in which case different sampled paths may start from different states or at
#' different times. By default, all sampled paths start from state 1 at time 0.
#' @param beta.state A matrix of dimension (no states) x (no transitions)
#' specifying estimated effects of times at which earlier states were reached
#' on subsequent transitions. If these are not in the model, the value
#' \code{NULL} (default) suffices; more information can be found under Details
#' @param clock Character argument, either \code{"forward"} (default) or
#' \code{"reset"}, specifying whether the time-scale of the cumulative hazards
#' is in forward time (\code{"forward"}) or duration in the present state
#' (\code{"reset"})
#' @param output One of \code{"state"}, \code{"path"}, or \code{"data"},
#' specifying whether states, paths, or data should be output.
#' @param tvec A numeric vector of time points at which the states or paths
#' should be evaluated. Ignored if \code{output}=\code{"data"}
#' @param cens An independent censoring distribution, given as a data frame
#' with time and Haz
#' @param M The number of sampled trajectories through the multi-state model.
#' The default is 10, since the procedure can become quite time-consuming
#' @param do.trace An integer, specifying that the replication number should be
#' written to the console every \code{do.trace} replications. Default is
#' \code{NULL} in which case no output is written to the console during the
#' simulation
#' @return M simulated paths through the multi-state model given by
#' \code{trans} and \code{Haz}. It is either a data frame with columns
#' \code{time}, \code{pstate1}, ..., \code{pstateS} for S states when
#' \code{output="state"}, or with columns \code{time}, \code{ppath1},...,
#' \code{ppathP} for the P paths specified in \code{\link{paths}}(trans) when
#' \code{output="path"}. When \code{output="data"}, the sampled paths are
#' stored in an \code{"msdata"} object, a data frame in long format such as
#' that obtained by \code{\link{msprep}}. This may be useful for
#' (semi-)parametric bootstrap procedures, in which case \code{cens} may be
#' used as censoring distribution (assumed to be independent of all transition
#' times and independent of any covariates).
#' @author Marta Fiocco, Hein Putter \email{H.Putter@@lumc.nl}
#' @references Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank
#' proportional hazards regression and simulation-based prediction for
#' multi-state models. \emph{Statistics in Medicine} \bold{27}, 4340--4358.
#' @keywords datagen
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # data in wide format, for transition 1 this is dataset E1 of
#' # Therneau & Grambsch (T&G)
#' tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'         dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'         x1=c(1,1,1,0,0,0),x2=c(6:1))
#' # data in long format using msprep
#' tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
#'         data=tg,keep=c("x1","x2"),trans=tmat)
#' # expanded covariates
#' tglong <- expand.covs(tglong,c("x1","x2"))
#' # Cox model with different covariate
#' cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
#'         data=tglong,method="breslow")
#' # new data, to check whether results are the same for transition 1 as T&G
#' newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
#' fit <- msfit(cx,newdata,trans=tmat)
#' tv <- unique(fit$Haz$time)
#' # mssample
#' set.seed(1234)
#' mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100)
#' set.seed(1234)
#' paths(tmat)
#' mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100,output="path")
#' set.seed(1234)
#' mssample(Haz=fit$Haz,trans=tmat,tvec=tv,M=100,output="data",do.trace=25)
#' 
#' @export mssample
`mssample` <- function(Haz, trans, history=list(state=1,time=0,tstate=NULL),
                    beta.state=NULL,
                    clock=c("forward","reset"),
                    output=c("state","path","data"), tvec, cens=NULL, M=10, do.trace=NULL)
{ 
    ###     Function to sample a path in the multi-state model
    ###     Input:
    ###         Haz: specific cumulative hazard for patient "i"
    ###         trans: matrix describing the multi-state model
    ###         history: the event-history of the patient, given by state and time
    ###         clock: either "forward" or "reset"
    ###         output: "state", "path", or "data"
    ###         tvec: vector of times
    ###         cens: distribution of censored times
    ###         M: number of replications
    ###      Output:
    ###         M simulated path in the muti-state model in the desired output format
    ###
    ### Set up
    output <- match.arg(output)
    clock <- match.arg(clock)
    K <- dim(trans)[1]
    trans2 <- to.trans2(trans)
    ntrans <- nrow(trans2)
    if (length(history$state)==1) history$state <- rep(history$state,M)
    if (length(history$time)==1) history$time <- rep(history$time,M)
    if (length(history$state)!=length(history$time)) stop("lengths of history$state and history$time differ")
    if (!is.null(history$tstate)) { # then should be either length K or dim 
        if (is.vector(history$tstate)) 
            if (length(history$tstate) != K) stop("length of history$tstate should equal no of states")
            else history$tstate <- matrix(history$tstate,K,M)
        if (is.null(beta.state)) stop("beta.state should be specified when history$tstate not null")
    }
    if (!is.null(beta.state))
        if (any(dim(beta.state) != c(K,ntrans)))
            stop("incorrect dimension of beta.state")
    if (output=="state") # to contain sum
        res <- matrix(0, length(tvec), K)
    else if (output=="path") {# to contain sum
        thepaths <- paths(trans)
        L <- nrow(thepaths)
        res <- matrix(0, length(tvec), L)
    }
    else res <- NULL
    for (m in 1:M) {
        if (!is.null(history$tstate))
            res1 <- mssample1(Haz, trans, history=list(state=history$state[m],
                        time=history$time[m], tstate=history$tstate[,m]), beta.state=beta.state, clock=clock, output=output,
                        tvec=tvec, cens=cens)
        else
            res1 <- mssample1(Haz, trans, history=list(state=history$state[m],
                        time=history$time[m], tstate=rep(0,K)), beta.state=beta.state, clock=clock, output=output,
                        tvec=tvec, cens=cens)
        if (output=="data") {
            res1[,1] <- m
            res <- rbind(res,res1)
        }
        else res <- res + res1
        if (!is.null(do.trace)) if (m %% do.trace == 0) {
            cat("Replication",m,"finished at",date(),"\n")
            flush.console()
        }
    }
    if (output=="state") {
        res <- data.frame(cbind(tvec,res/M))
        names(res) <- c("time",paste("pstate",1:K,sep=""))
    }
    else if (output=="path") {
        res <- data.frame(cbind(tvec,res/M))
        names(res) <- c("time",paste("ppath",1:L,sep=""))
    }
    else if (output=="data") {
        res <- data.frame(res)
        names(res) <- c("id","Tstart","Tstop","duration","from","to","status","trans")
        attr(res, "trans") <- trans
        class(res) <- "msdata"
    }
    return(res)  
}

`mssample1` <- function(Haz, trans, history, beta.state, clock, output, tvec, cens)
{ 
###
###     Function to sample a single path in the multi-state model
###     Used internally in function mssample and msboot
###
    ### First sample (once) from censoring distribution, this will be the follow-up time 
    if (!is.null(cens)) {
        pcens <- diff(c(0,1-cens$surv))
        idx <- sample(1:length(cens$time), size=1, prob=pcens)
        fut <- cens$time[idx]
        censtime <- list(time=fut, jump=ifelse(idx>1,cens$Haz[idx]-cens$Haz[idx-1],cens$Haz[idx]))
    }
    else censtime <- NULL
    K <- dim(trans)[1]
    trans2 <- to.trans2(trans)
    from <- to <- history$state
    tcond <- t0 <- Tstart <- history$time
    if (output=="state")
        res <- matrix(0, length(tvec), K)
    else if (output=="path") {
        thepaths <- paths(trans)
        path <- c(to, rep(NA,ncol(thepaths)-1))
        res <- matrix(0, length(tvec), nrow(thepaths))
   }
    else res <- NULL
    ### keep track of when states were visited
    tstates <- history$tstate
    while (!is.na(to)) {
        from <- to
        nstates <- trans[from,]
        transs <- nstates[!is.na(nstates)]
        allto <- which(!is.na(nstates))
        ntr <- length(transs)
              
        if (ntr!=0 ) { # if not yet in absorbing state
            transnos <- transs
            for (tr in 1:ntr)
                Haz$Haz[Haz$trans==transnos[tr]] <-
                  exp(sum(beta.state[,transnos[tr]]*tstates)) *
                    Haz$Haz[Haz$trans==transnos[tr]]
            whh <- which(!is.na(match(Haz$trans,transnos)))
            if (clock=="forward") {
                crs <- crsample(Haz[whh,], tcond, censtime)
                tcond <- Tstop <- crs$t
            }
            else {
                crs <- crsample(Haz[whh,], t0, censtime)
                t0 <- 0
                tcond <- Tstop <- crs$t + tcond
            }
            transno <- crs$trans
            if (is.na(transno)) to <- NA
            else {
                to <- trans2$to[transno]
                tstates[to] <- Tstop
            }
            if (output=="state") {
                res[((tvec>=Tstart)&(tvec<Tstop)),from] <- 1
                Tstart <- Tstop
            }
            else if (output=="path") {
                idx <- which(apply(thepaths,1,function(x) identical(x,path)))
                res[((tvec>=Tstart)&(tvec<Tstop)),idx] <- 1
                path[which(is.na(path))[1]] <- to
                Tstart <- Tstop
            }
            else {
                res1 <- matrix(c(rep(NA,ntr),rep(Tstart,ntr),rep(Tstop,ntr),rep(Tstop-Tstart,ntr),rep(from,ntr),allto,rep(0,2*ntr)),ntr,8)
                res1[res1[,6]==to,7] <- 1
                res1[,8] <- trans[from,allto] # trans
                Tstart <- Tstop
                res <- rbind(res,res1)
            }
        } 
        else { # in absorbing state:
            to <- NA
            if (output=="state") {
                res[tvec>=Tstart,from] <- 1
            }
            else if (output=="path") {
                idx <- which(apply(thepaths,1,function(x) identical(x,path)))
                res[tvec>=Tstart,idx] <- 1
                path[which(is.na(path))[1]] <- to
            }
            else {
                res1 <- matrix(c(rep(NA,ntr),rep(Tstart,ntr),rep(Tstop,ntr),rep(Tstop-Tstart,ntr),rep(from,ntr),allto,rep(0,2*ntr)),ntr,8)
                res1[res1[,6]==to,7] <- 1
                res1[,8] <- trans[from,allto] # trans
                res <- rbind(res,res1)
            }
        }
    }
    return(res)  
}

`crsample` <- function(Haz, tcond=0, censtime=NULL) ### censtime here is list with time and last jump (=pp.fut)
{
    if (is.null(censtime)) fut <- Inf
    else
        fut <- censtime$time
    transs <- Haz$trans
    transun <- unique(transs)
    K <- length(transun)
    tt <- sort(unique(Haz$time))
    n <- length(tt)
    cim <- matrix(NA, n, 3*K+4)
    ci <- as.data.frame(cim)
    names(ci)[1] <- "time"
    names(ci)[2:(K+1)] <- paste("Haz",as.character(1:K), sep="")
    names(ci)[(K+2):(2*K+1)] <- paste("haz",as.character(1:K), sep="")
    names(ci)[(2*K+2):(3*K+1)] <- paste("CI",as.character(1:K), sep="")
    names(ci)[3*K+2] <- "hazsum"
    names(ci)[3*K+3] <- "Hazsum"
    names(ci)[3*K+4] <- "S0"
    ci$time <- tt
    for (k in 1:K) { # select the elements for each transition in the block
        wh <- which(Haz$trans==transun[k])
        idx <- match(Haz$time[wh],tt)
        ci[,k+1][idx] <- Haz$Haz[wh]
        ci[,k+1] <- NAfix(ci[,k+1],subst=0)
        ci[,K+1+k] <- diff(c(0,ci[,k+1])) # the hazard for transition k
    }
    ## select only t > tcond, and adjust the cumulative hazards
    ci <- ci[ci$time>tcond,]
    n <- nrow(ci)
    for (k in 1:K) # each cumulative hazard is adjusted (cumulative sum of hazard)
        ci[,k+1] <- cumsum(ci[,K+1+k])
    ### compute baseline survival S0, need to sum all columns Hazards
    if (K==1) ci$hazsum <- ci[,3]
    else ci$hazsum <- apply(ci[,((K+2):(2*K+1))],1,sum)
    ci$S0 <- cumprod(1-ci$hazsum)
    ci$Hazsum <- -log(ci$S0)
    nci <- nrow(ci)
    ### Follow Dabrowska in first sampling from the sum of the cumulative hazards
    k <- NA
    tsample <- Hazsample(data.frame(time=ci$time, Haz=ci$Hazsum))
    ### If follow-up time is less than tsample, no transition is made
    if (fut < tsample) crt <- fut
    else { # else decide between transitions
        crt <- tsample
        if (fut>tsample) # the size of the hazard jumps of the transitions
                         # determine the probabilities of choosing the transitions
        {
            k <- sample(1:K,size=1,prob=ci[which(ci$time==tsample),(K+2):(2*K+1)])
        }
        else # i.e. if there is fut==tsample, then the jump in the censoring distribution
             # also plays in the lottery, unless both are Inf, then it doesn't matter
            if (crt!=Inf) {
                k <- sample(c(1:K,NA),size=1,
                        prob=c(ci[which(ci$time==tsample),(K+2):(2*K+1)],censtime$jump))
            }
    }
    if (!is.na(k)) trans <- unique(Haz$trans)[k] # return transition k
    else trans <- NA
    return(list(t=crt,trans=trans))
}

`Hazsample` <- function(Haz, size=1, replace=TRUE)
{
    ### Function to sample from a distribution given by a cumulative hazard
    ### Input:
    ###     time: vector containing time points to sample from
    ###     H: the cumulative hazard at those time points
    ###     size: the size of the sample
    ###     replace: sample with or without replacement? (default=TRUE, with replacement)
    ### Output:
    ###     A random realisation from H
    p <- diff(c(0,1-exp(-Haz$Haz)))
    p <- c(p, exp(-Haz$Haz[nrow(Haz)])) # add probability of sampling time=Inf
    return(sample(c(Haz$time, Inf), size=size, prob=p, replace=replace))
}
