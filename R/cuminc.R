`Cuminc` <- function(time, status, data, group, failcodes, na.status=c("remove","extra"), variance=TRUE)
{
    ## time
    if (!is.vector(time)) stop("argument \"time\" not of correct type")
    if (is.character(time)) {
        if (length(time) != 1)
            stop("single character string required for \"time\" argument")
        if (missing(data)) stop("argument \"data\" missing")
        time <- data[[time]]
    }
    ## status
    if (!is.vector(status)) stop("argument \"status\" not of correct type")
    if (is.character(status)) {
        if (length(status) == 1) {
            if (missing(data)) stop("argument \"data\" missing")
            status <- data[[status]]
        }
    }
    ## check for compatibility of time and status
    if (length(status) != length(time)) stop("lengths of time and status do not match")
    na.status <- match.arg(na.status)
    if (na.status=="remove") wh <- which(!is.na(status))
    time <- time[wh]; status <- status[wh]
    n <- length(time)
    ## For future calculations it is easiest to work with status variable
    ## with values 1 through K for K causes of failure and 0 for censorings
    if (missing(failcodes)) {
        failcodes <- sort(unique(status),na.last=TRUE)
        if (!is.factor(status)) {
            failcodes <- sort(unique(status),na.last=TRUE)
            if (failcodes[1] < 0) stop("smallest value should be 0 to indicate censoring")
            else if (failcodes[1]==0) failcodes <- failcodes[-1]
        }
        K <- length(failcodes)
        status <- match(status,failcodes)
        status[is.na(status)] <- 0
    }
    else {
        status <- match(status,failcodes)
        status[is.na(status)] <- 0
        tbl <- table(status)
        fc <- as.numeric(names(tbl))
        if (fc[1]==0) fc <- fc[-1]
        failcodes <- failcodes[fc]
        K <- length(failcodes)
        status[status>0] <- match(status[status>0],fc)
    }
    if (K==0) stop("Only censored observations found")
    if (K==1) stop("Only 1 failure cause found; please use survfit instead")
    if (any(time<0)) stop("Negative time values not allowed")
    ## Make transition matrix and prepare data in long format
    tmat <- trans.comprisk(K)
    ttime <- matrix(time,n,K+1) # time also in first row, doesn't matter
    sstatus <- matrix(0,n,K+1)
    for (k in 1:K) sstatus[status==k,k+1] <- 1
    if (missing(group)) {
        msdata <- msprep(time=ttime,status=sstatus,trans=tmat)
        c0 <- coxph(Surv(time,status)~strata(trans), data=msdata, method="breslow")
        msf0 <- msfit(c0,variance=variance,vartype="greenwood",trans=tmat)
        if (any(msf0$Haz$Haz[msf0$Haz$time==0]>0)) { # make sure that events at t=0 are correctly dealt with
            ptg <- probtrans(msf0,predt=-0.0001,variance=variance,method="greenwood")[[1]]
            ptg$time[1] <- 0
        }
        else ptg <- probtrans(msf0,predt=0,variance=variance,method="greenwood")[[1]]
        names(ptg)[2:(K+2)] <- c("Surv",paste("CI",failcodes,sep="."))
        if (variance) names(ptg)[(K+3):(2*K+3)] <- paste("se",names(ptg)[2:(K+2)],sep="")
    }
    else {
        ## group (bit of code taken and adapted from msprep, perhaps not all necessary)
        if (!(is.matrix(group)|(is.data.frame(group)))) { # then group should be a vector
            if (is.character(group)) { # character vector
                if (missing(data)) stop("argument \"data\" is missing, with no default")
                ngroup <- length(group)
                if (ngroup!=1) stop("only single grouping variable possible")
                kcols <- match(group,names(data))
                if (any(is.na(kcols))) stop("\"group\" column not in data")
                groupname <- group
                group <- data[,kcols]
            }
            else { # other vector, should be of length n
                ngroup <- 1
                groupname <- names(group)
                if (length(group) != n) stop("argument \"group\" has incorrect dimension")
            }
        }
        else {
            ngroup <- ncol(group)
            groupname <- names(group)
            if (nrow(group) != n) stop("argument \"group\" has incorrect dimension")
            if (ngroup!=1) stop("only single grouping variable possible")
            group <- group[,1] # coerce to vector
        }
        if (is.null(groupname)) groupname <- "group"
        ## msprep and the rest
        msdata <- msprep(time=ttime,status=sstatus,trans=tmat,keep=group)
        names(msdata)[ncol(msdata)] <- "group"
        groupvals <- sort(unique(msdata$group))
        G <- length(groupvals)
        ptg <- NULL
        for (g in 1:G) {
            msdatag <- msdata[msdata$group==groupvals[g],]
            tbl <- table(msdatag$status,msdatag$trans)
            if (dim(tbl)[1]>1) # not if there are only censored observations for this group
            {
                ### If for a certain group value one or more of the causes of failure is not present
                ### probtrans will work, but pstates are assigned to 1 to number of valid causes
                ### then these valid causes should be reassigned to the correct causes
                ### In case of only one cause of failure with events, call survfit instead
                wh <- which(tbl[2,]>0) # these should be the causes that are present
                nwh <- length(wh)
                msdatag <- msdatag[msdatag$trans %in% wh,]
                tbl <- table(msdatag$status,msdatag$trans)
                if (dim(tbl)[2] == 1) { # only one cause of failure with events
                    fit <- survfit(Surv(time, status) ~ 1, data=msdatag)
                    ptgg <- data.frame(time=fit$time,surv=fit$surv,CI=1-fit$surv,seSurv=fit$std.err)
                    ptgg$seSurv <- ptgg$surv * ptgg$seSurv # survfit gives SE of cumulative hazard, want that of survival
                    ptgg$seSurv[ptgg$surv==0] <- 0
                    ptgg$seCI <- ptgg$seSurv
                }
                else {                
                    c0 <- coxph(Surv(time,status)~strata(trans), data=msdatag, method="breslow")
                    msf0 <- msfit(c0,variance=variance,vartype="greenwood",trans=tmat)
                    if (any(msf0$Haz$Haz[msf0$Haz$time==0]>0)) {
                        ptgg <- probtrans(msf0,predt=-0.0001,variance=variance,method="greenwood")[[1]]
                        ptgg$time[1] <- 0
                    }
                    else ptgg <- probtrans(msf0,predt=0,variance=variance,method="greenwood")[[1]]
                }
                if (nwh==1) { # only one cause of failure with events
                    tmp <- ptgg
                    if (variance) ptgg <- matrix(NA,nrow(tmp),2*K+3)
                    else ptgg <- matrix(NA,nrow(tmp),K+3)
                    ptgg <- as.data.frame(ptgg)
                    ptgg[,1:2] <- tmp[,1:2]
                    ptgg[,2 + wh] <- tmp[,2 + 1:length(wh)]
                    ptgg[,(3:(K+2))[-wh]] <- 0
                    if (variance) {
                        ptgg[,K+3] <- tmp[,nwh+3]
                        ptgg[,K+3 + wh] <- tmp[,nwh+3 + 1:length(wh)]
                        ptgg[,((K+4):(2*K+3))[-wh]] <- 0
                    }
                }
                else if (nwh < K) { # more than one, fewer than K causes of failure with events
                    ptgg[,2 + wh] <- ptgg[,2 + 1:nwh]
                    ptgg[,(3:(K+2))[-wh]] <- 0
                    if (variance) {
                        ptgg[,K+3 + wh] <- ptgg[,K+3 + 1:nwh]
                        ptgg[,((K+4):(2*K+3))[-wh]] <- 0
                    }
                }
                names(ptgg)[1:(K+2)] <- c("time","Surv",paste("CI",failcodes,sep="."))
                if (variance) names(ptgg)[(K+3):(2*K+3)] <- paste("se",names(ptgg)[2:(K+2)],sep="")
                ptgg$group <- g
                ptg <- rbind(ptg,ptgg)
            }
        }
        tbl <- table(group)
        if (is.factor(group)) {
            wh <- which(tbl>0) # the causes that are present
            ptg$group <- factor(ptg$group,levels=1:G,labels=levels(group)[wh])
        } else {
            wh <- as.numeric(names(tbl))
            ptg$group <- wh[ptg$group]
        }
        # Put group first and rename to original variable
        nc <- ncol(ptg)
        ptg <- ptg[,c(nc,1:(nc-1))]
        names(ptg)[1] <- groupname
    }
    return(ptg)
}
