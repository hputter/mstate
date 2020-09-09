#' Compute subject-specific or overall transition probabilities with standard
#' errors
#' 
#' This function computes subject-specific or overall transition probabilities
#' in multi-state models. If requested, also standard errors are calculated.
#' 
#' For details refer to de Wreede, Fiocco & Putter (2010).
#' 
#' @param object \link{msfit} object containing estimated cumulative hazards
#' for each of the transitions in the multi-state model and, if standard errors
#' are requested, (co)variances of these cumulative hazards for each pair of
#' transitions
#' @param predt A positive number indicating the prediction time. This is
#' either the time at which the prediction is made (if \code{direction}=
#' \code{"forward"}) or the time for which the prediction is to be made (if
#' \code{direction}=\code{"fixedhorizon"})
#' @param direction One of \code{"forward"} (default) or \code{"fixedhorizon"},
#' indicating whether prediction is forward or for a fixed horizon
#' @param method A character string specifying the type of variances to be
#' computed (so only needed if either \code{variance} or \code{covariance} is
#' \code{TRUE}). Possible values are \code{"aalen"} or \code{"greenwood"}
#' @param variance Logical value indicating whether standard errors are to be
#' calculated (default is \code{TRUE})
#' @param covariance Logical value indicating whether covariances of transition
#' probabilities for different states are to be calculated (default is
#' \code{FALSE})
#' @return An object of class \code{"probtrans"}, which is a list of which item
#' [[s]] contains a data frame with the estimated transition probabilities (and
#' standard errors if \code{variance}=\code{TRUE}) from state s.  If
#' \code{covariance}=\code{TRUE}, item \code{varMatrix} contains an array of
#' dimension K^2 x K^2 x (nt+1) (with K the number of states and nt the
#' distinct transition time points); the time points correspond to those in the
#' data frames with the estimated transition probabilities.  Finally, there are
#' items \code{trans}, \code{method}, \code{predt}, \code{direction}, recording
#' the transition matrix, and the method, predt and direction arguments used in
#' the call to probtrans.  Plot and summary methods have been defined for
#' \code{"probtrans"} objects.
#' @author Liesbeth de Wreede and Hein Putter \email{H.Putter@@lumc.nl}
#' @references Andersen PK, Borgan O, Gill RD, Keiding N (1993).
#' \emph{Statistical Models Based on Counting Processes}. Springer, New York.
#' 
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: Competing
#' risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
#' 2389--2430.
#' 
#' Therneau TM, Grambsch PM (2000). \emph{Modeling Survival Data: Extending the
#' Cox Model}. Springer, New York.
#' 
#' de Wreede LC, Fiocco M, and Putter H (2010). The mstate package for
#' estimation and prediction in non- and semi-parametric multi-state and
#' competing risks models. \emph{Computer Methods and Programs in Biomedicine}
#' \bold{99}, 261--274.
#' 
#' de Wreede LC, Fiocco M, and Putter H (2011). mstate: An R Package for the
#' Analysis of Competing Risks and Multi-State Models. \emph{Journal of
#' Statistical Software}, Volume 38, Issue 7.
#' @keywords survival
#' @examples
#' 
#' # transition matrix for illness-death model
#' tmat <- trans.illdeath()
#' # data in wide format, for transition 1 this is dataset E1 of
#' # Therneau & Grambsch (2000)
#' tg <- data.frame(illt=c(1,1,6,6,8,9),ills=c(1,0,1,1,0,1),
#'         dt=c(5,1,9,7,8,12),ds=c(1,1,1,1,1,1),
#'         x1=c(1,1,1,0,0,0),x2=c(6:1))
#' # data in long format using msprep
#' tglong <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
#' 		data=tg,keep=c("x1","x2"),trans=tmat)
#' # events
#' events(tglong)
#' table(tglong$status,tglong$to,tglong$from)
#' # expanded covariates
#' tglong <- expand.covs(tglong,c("x1","x2"))
#' # Cox model with different covariate
#' cx <- coxph(Surv(Tstart,Tstop,status)~x1.1+x2.2+strata(trans),
#'         data=tglong,method="breslow")
#' summary(cx)
#' # new data, to check whether results are the same for transition 1 as
#' # those in appendix E.1 of Therneau & Grambsch (2000)
#' newdata <- data.frame(trans=1:3,x1.1=c(0,0,0),x2.2=c(0,1,0),strata=1:3)
#' HvH <- msfit(cx,newdata,trans=tmat)
#' # probtrans
#' pt <- probtrans(HvH,predt=0)
#' # predictions from state 1
#' pt[[1]]
#' 
#' @export probtrans
`probtrans` <- function(object,predt,direction=c("forward","fixedhorizon"),
                method=c("aalen","greenwood"),variance=TRUE,covariance=FALSE)
{
    if (!inherits(object, "msfit"))
        stop("'object' must be a 'msfit' object")
    method <- match.arg(method)
    direction <- match.arg(direction)
    trans <- object$trans
    transit <- to.trans2(trans)
    numtrans <- nrow(transit)
    stackhaz <- object$Haz
    stackvarhaz <- object$varHaz
    for (i in 1:numtrans)
        stackhaz$dhaz[stackhaz$trans==i] <- diff(c(0,stackhaz$Haz[stackhaz$trans==i]))
    if (direction=="forward")
        stackhaz <- stackhaz[stackhaz$time > predt,]
    else stackhaz <- stackhaz[stackhaz$time <= predt,]

    untimes <- sort(unique(stackhaz$time))
    TT <- length(untimes)
    S <- nrow(trans)

    if (covariance) variance <- TRUE # if covariance=TRUE and variance=FALSE, variance is overruled
    if (direction=="forward") {
        if (variance==TRUE) res <- array(0,c(TT+1,2*S+1,S)) # 2*S+1 for time, probs (S), se (S)
        else res <- array(0,c(TT+1,S+1,S)) # S+1 for time, probs (S)
        # first line (in case of forward) contains values for t=predt
        res[1,1,] <- predt
        for (j in 1:S) res[1, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
        if (variance) res[1,(S+2):(2*S+1),] <- 0
    }
    else {
        # situation for backward is different from forward,
        # depends on whether predt is an event time
        if (predt %in% untimes) {
            if (variance) res <- array(0,c(TT+1,2*S+1,S)) # t=0 and all event times
            else res <- array(0,c(TT+1,S+1,S))
            res[TT+1,1,] <- predt
            for (j in 1:S) res[TT+1, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
            if (variance) res[TT+1,(S+2):(2*S+1),] <- 0
        }
        else {
            if (variance) res <- array(0,c(TT+2,2*S+1,S)) # t=0, event times, t=predt
            else res <- array(0,c(TT+2,S+1,S))
            res[TT+1,1,] <- max(untimes)
            for (j in 1:S) res[TT+1, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
            if (variance) res[TT+1,(S+2):(2*S+1),] <- 0
            res[TT+2,1,] <- predt
            for (j in 1:S) res[TT+2, 1+j,] <- rep(c(0,1,0), c(j-1,1,S-j))
            if (variance) res[TT+2,(S+2):(2*S+1),] <- 0
        }
    }
    P <- diag(S)
    
    if (covariance) {
      varParr <- array(0,c(S^2,S^2,TT+1))
      if ((direction=="fixedhorizon") & !(predt %in% untimes))
        varParr <- array(0,c(S^2,S^2,TT+2))
      ffrom <- rep(1:S, S)
      tto <- rep(1:S, rep(S,S))
      fromto <- paste("from",ffrom,"to",tto,sep="")
      if (direction=="forward")
        dimnames(varParr) <- list(fromto,fromto,c(predt,untimes))
      else {
        if (predt %in% untimes) dimnames(varParr) <- list(fromto,fromto,c(0,untimes))
        else dimnames(varParr) <- list(fromto,fromto,c(0,untimes,predt))
      }
    }
    if (variance) {
        varP <- matrix(0,S^2,S^2)
        if (direction=="forward") {
            varAnew <- array(0,c(S,S,S,S))
            if (predt !=0) {
                tmin <- max(stackvarhaz$time[stackvarhaz$time <=predt])
                varHaz <- stackvarhaz[stackvarhaz$time==tmin,]
                lHaz <- nrow(varHaz)
                for (j in 1:lHaz)
                {
                    from1 <- transit$from[transit$transno==varHaz$trans1[j]]
                    to1 <- transit$to[transit$transno==varHaz$trans1[j]]
                    from2 <- transit$from[transit$transno==varHaz$trans2[j]]
                    to2 <- transit$to[transit$transno==varHaz$trans2[j]]
                    varAnew[from1, to1, from2, to2] <-
                        varAnew[from2, to2, from1, to1] <- varHaz$varHaz[j]
                }
            }
        }
        else { 
            # varA on last event time (only those elements borrowed
            # directly from object, rest is calculated later)
            varA <- array(0,c(S,S,S,S))
            varHaz <- stackvarhaz[stackvarhaz$time==untimes[TT],]
            lHaz <- nrow(varHaz)
            for (j in 1:lHaz)
            {
                from1 <- transit$from[transit$transno==varHaz$trans1[j]]
                to1 <- transit$to[transit$transno==varHaz$trans1[j]]
                from2 <- transit$from[transit$transno==varHaz$trans2[j]]
                to2 <- transit$to[transit$transno==varHaz$trans2[j]]
                varA[from1, to1, from2, to2] <-
                    varA[from2, to2, from1, to1] <- varHaz$varHaz[j]
            }
        }
    }
    
    for (i in 1:TT)
    {
        idx <- ifelse(direction=="forward",i,TT+1-i)
        tt <- untimes[idx]
        Haztt <- stackhaz[stackhaz$time==tt,]
        lHaztt <- nrow(Haztt)
        # build S x S matrix IplusdA
        IplusdA <- diag(S)
        for (j in 1:lHaztt)
        {
            from <- transit$from[transit$transno==Haztt$trans[j]]
            to <- transit$to[transit$transno==Haztt$trans[j]]
            IplusdA[from, to] <- Haztt$dhaz[j]
            IplusdA[from, from] <- IplusdA[from, from] - Haztt$dhaz[j]
        }
        if (any(diag(IplusdA)<0))
            warning("Warning! Negative diagonal elements of (I+dA); the estimate may not be meaningful. \n")
        
        if (variance) {
            if (direction=="forward"){
                varA <- varAnew
                varAnew <- array(0,c(S,S,S,S))
                varHaztt <- stackvarhaz[stackvarhaz$time==tt,]
                lHaztt <- nrow(varHaztt)
                for (j in 1:lHaztt)
                {
                    from1 <- transit$from[transit$transno==varHaztt$trans1[j]]
                    to1 <- transit$to[transit$transno==varHaztt$trans1[j]]
                    from2 <- transit$from[transit$transno==varHaztt$trans2[j]]
                    to2 <- transit$to[transit$transno==varHaztt$trans2[j]]
                    varAnew[from1, to1, from2, to2] <-
                        varAnew[from2, to2, from1, to1] <- varHaztt$varHaz[j]
                }
                vardA <- varAnew - varA
            }
            else {
                varAttmin <- array(0,c(S,S,S,S))
                varHazttmin <- stackvarhaz[stackvarhaz$time==untimes[idx-1],]
                lHazttmin <- nrow(varHazttmin)
                for (j in 1:lHazttmin)
                {
                    from1 <- transit$from[transit$transno==varHazttmin$trans1[j]]
                    to1 <- transit$to[transit$transno==varHazttmin$trans1[j]]
                    from2 <- transit$from[transit$transno==varHazttmin$trans2[j]]
                    to2 <- transit$to[transit$transno==varHazttmin$trans2[j]]
                    varAttmin[from1, to1, from2, to2] <-
                        varAttmin[from2, to2, from1, to1] <- varHazttmin$varHaz[j]
                }
                vardA <- varA - varAttmin
                varA <- varAttmin # ready for the next round
            }
        
            for (from in 1:S)
            {
                for (from2 in 1:S)
                {
                    for (to2 in 1:S)
                    {
                        if (to2!=from2)
                            vardA[from, from, from2, to2] <-
                                vardA[from2, to2, from, from] <-
                                    -sum(vardA[from,-from,from2,to2])
                    }
                }
            }
            for (from in 1:S)
            {
                for (from2 in 1:S)
                    vardA[from,from,from2,from2] <-
                        vardA[from2,from2,from,from] <-
                        -sum(vardA[from,from,from2,-from2])
            }
            vardA <- matrix(vardA,S^2,S^2)
        }
        
        if (method=="aalen")
        {
            if (direction=="forward") 
            {  
                P <- P %*% IplusdA
                if (variance) {
                    tmp1 <- kronecker(t(IplusdA),diag(S)) %*% varP %*% kronecker(IplusdA,diag(S))
                    tmp2 <- kronecker(diag(S),P) %*% vardA %*% kronecker(diag(S),t(P))
                    varP <- tmp1 + tmp2
                }
            }
            else {
                if (variance) {
                    tmp1 <- kronecker(diag(S), IplusdA) %*% varP %*% kronecker(diag(S), t(IplusdA))
                    tmp2 <- kronecker(t(P), IplusdA) %*% vardA %*% kronecker(P,t(IplusdA))
                    varP <- tmp1+tmp2
                }
                P <- IplusdA %*% P
            }
        }
        if (method=="greenwood")
        {
            if (direction=="forward")
            {
                if (variance) {
                    tmp1 <- kronecker(t(IplusdA),diag(S)) %*% varP %*% kronecker(IplusdA,diag(S))
                    tmp2 <- kronecker(diag(S),P) %*% vardA %*% kronecker(diag(S),t(P))
                    varP <- tmp1 + tmp2
                }
                P <- P %*% IplusdA
            }
            else {
                if (variance) {
                    tmp1 <- kronecker(diag(S), IplusdA) %*% varP %*% kronecker(diag(S), t(IplusdA))
                    tmp2 <- kronecker(t(P), diag(S)) %*% vardA %*% kronecker(P, diag(S))
                    varP <- tmp1+tmp2
                }
                P <- IplusdA %*% P
            }
        }
        if (variance) {
            seP <- sqrt(diag(varP))
            seP <- matrix(seP,S,S)
        }
        if (covariance) {
          if (direction=="forward") varParr[,,i+1] <- varP
          else {
            varParr[,,idx+1] <- varP
          }
        }
        if (direction=="forward") {
            res[idx+1,1,] <- tt
            res[idx+1,2:(S+1),] <- t(P)
            if (variance) res[idx+1,(S+2):(2*S+1),] <- t(seP)
        }
        else {
            res[idx,1,] <- ifelse(i==TT,0,untimes[TT-i])
            res[idx,2:(S+1),] <- t(P)
            if (variance) res[idx,(S+2):(2*S+1),] <- t(seP)
        }
    }
    if (covariance & (direction=="fixedhorizon"))
      varParr[,,1] <- varParr[,,2]
    ### res[,,s] contains prediction from state s, convert to list of dataframes
    res2  <- vector("list", S)
    for (s in 1:S) {
        tmp <- as.data.frame(res[,,s])
        if (min(dim(tmp))==1) tmp <- res[,,s]
        if (variance) names(tmp) <- c("time",paste("pstate",1:S,sep=""),paste("se",1:S,sep=""))
        else names(tmp) <- c("time",paste("pstate",1:S,sep=""))
        res2[[s]] <- tmp
    }
    if (covariance) res2$varMatrix <- varParr
    res2$trans <- trans
    res2$method <- method
    res2$predt <- predt
    res2$direction <- direction
    class(res2) <- "probtrans"
    return(res2)
}
