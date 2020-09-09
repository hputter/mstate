`redrank.iter` <- function(data, cols, colsR, fullcols, R,
                        clock, strata, Gamma.iter, print.level = 0)
{
    if (!inherits(data, "msdata"))
        stop("'data' must be an 'msdata' object")
    trans <- attr(data, "trans")
    data[,c(cols,colsR,fullcols)] <- unlist(data[,c(cols,colsR,fullcols)])
    nfull <- length(fullcols)
    covariates <- names(data)[cols]
    fullcovariates <- names(data)[fullcols]
    n <- nrow(data)
    K <- max(trans,na.rm=TRUE)
    p <- length(cols)
    for (k in 1:K) { ### W[k,r] = Gamma.iter[r,k] * Z
        wh <- which(data$trans == k)
        for (r in 1:R) {
            data[wh, colsR[r,]] <- Gamma.iter[r,k] * data[wh, colsR[r,]]
        }
    }
    data$trans <- factor(data$trans) ### lost in the unlist above
    covs.R <- names(data)[as.vector(t(colsR))]
    if (is.null(strata)) strata <- "trans"
    strata.expr <- paste("strata(",strata,")+")
    if (clock=="forward")
        expr1 <- paste("cox.itr1 <- coxph(Surv(Tstart,Tstop,status)~",
            strata.expr,
            paste(covs.R, collapse = "+"),sep="")
    else
        expr1 <- paste("cox.itr1 <- coxph(Surv(time,status)~",
            strata.expr,
            paste(covs.R, collapse = "+"),sep="")
    if (!is.null(fullcols)) expr1 <- paste(expr1,
        "+", paste(fullcovariates, collapse = "+"),sep="")
    expr1 <- paste(expr1,", data=data, na.action=\"na.exclude\")", sep = "")
    cox.itr1 <- eval(parse(text = expr1, n = 1))
    if (print.level > 0) {
        cat("\nStratified Cox regression:\n\n")
        print(cox.itr1)
    }
    loglik <- cox.itr1$loglik[2]
    if (print.level > -1) cat("\nAlpha =", round(cox.itr1$coef, 5))
    Alpha <- matrix(cox.itr1$coef[1:(p*R)],p,R)
    Beta2 <- cox.itr1$coef[-(1:(p*R))]
    ncd <- ncol(data)
    data <- cbind(data,as.matrix(data[, cols]) %*% Alpha)
    AlphaX.R <- paste("AlphaX",as.character(1:R),sep="")
    names(data)[((ncd+1):(ncd+R))] <- AlphaX.R
    attr(data, "trans") <- trans
    class(data) <- c("msdata", "data.frame")
    data <- expand.covs(data,AlphaX.R)
    AlphaX.RK <- names(data)[((ncd+R+1):(ncd+R+R*K))] 
    if (clock=="forward")
        expr2 <- paste("cox.itr2 <- coxph(Surv(Tstart,Tstop,status)~",
            strata.expr,
            paste(AlphaX.RK, collapse = "+"),sep="")
    else
        expr2 <- paste("cox.itr2 <- coxph(Surv(time,status)~",
            strata.expr,
            paste(AlphaX.RK, collapse = "+"),sep="")
    if (!is.null(fullcols)) expr2 <- paste(expr2,
        "+", paste(fullcovariates, collapse = "+"),sep="")
    expr2 <- paste(expr2,", data=data, na.action=\"na.exclude\")", sep = "")
    cox.itr2 <- eval(parse(text = expr2, n = 1))
    if (print.level > 0) {
        cat("\n\nCox regression on scores\n\n")
        print(cox.itr2)
    }   
    Gamma.iter <- t(matrix(cox.itr2$coef[1:(K*R)],K,R))
    Beta2 <- cox.itr2$coef[-(1:(K*R))]
    return(list(Gamma = Gamma.iter, Alpha = Alpha, Beta2 = Beta2,
                loglik = loglik, cox.itr1 = cox.itr1))
}





#' Reduced rank proportional hazards model for competing risks and multi-state
#' models
#' 
#' This function estimates regression coefficients in reduced rank proportional
#' hazards models for competing risks and multi-state models.
#' 
#' For details refer to Fiocco, Putter & van Houwelingen (2005, 2008).
#' 
#' @param redrank Survival formula, starting with either Surv(time,status) ~ or
#' with Surv(Tstart,Tstop,status) ~, followed by a formula containing
#' covariates for which a reduced rank estimate is to be found
#' @param full Optional, formula specifying that part which needs to be
#' retained in the model (so not subject to reduced rank)
#' @param data Object of class 'msdata', as prepared for instance by
#' \code{\link{msprep}}, in which to interpret the \code{redrank} and,
#' optionally, the \code{full} formulas
#' @param R Numeric, indicating the rank of the solution
#' @param strata Name of covariate to be used inside the
#' \code{\link[survival:strata]{strata}} part of
#' \code{\link[survival:coxph]{coxph}}
#' @param Gamma.start A matrix of dimension K x R, with K the number of
#' transitions and R the rank, to be used as starting value
#' @param method The method for handling ties in
#' \code{\link[survival:coxph]{coxph}}
#' @param eps Numeric value determining when the iterative algorithm is
#' finished (when for two subsequent iterations the difference in
#' log-likelihood is smaller than \code{eps})
#' @param print.level Determines how much output is written to the screen; 0:
#' no output, 1: iterations, for each iteration solutions of Alpha, Gamma,
#' log-likelihood, 2: also the Cox regression results
#' @return A list with elements \item{Alpha}{the Alpha matrix} \item{Gamma}{the
#' Gamma matrix} \item{Beta}{the Beta matrix corresponding to
#' \code{covariates}} \item{Beta2}{the Beta matrix corresponding to
#' \code{fullcovs}} \item{cox.itr1}{the \code{\link[survival:coxph]{coxph}}
#' object resulting from the last call giving \code{Alpha}} \item{alphaX}{the
#' matrix of prognostic scores given by \code{Alpha}, n x R, with n number of
#' subjects} \item{niter}{the number of iterations needed to reach convergence}
#' \item{df}{the number of regression parameters estimated} \item{loglik}{the
#' log-likelihood}
#' @author Marta Fiocco and Hein Putter \email{H.Putter@@lumc.nl}
#' @references Fiocco M, Putter H, van Houwelingen JC (2005). Reduced rank
#' proportional hazards model for competing risks. \emph{Biostatistics}
#' \bold{6}, 465--478.
#' 
#' Fiocco M, Putter H, van Houwelingen HC (2008). Reduced-rank proportional
#' hazards regression and simulation-based prediction for multi-state models.
#' \emph{Statistics in Medicine} \bold{27}, 4340--4358.
#' 
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: Competing
#' risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
#' 2389--2430.
#' @keywords survival
#' @examples
#' 
#' \dontrun{
#'     # This reproduces the results in Fiocco, Putter & van Houwelingen (2005)
#'     # Takes a while to run
#'     data(ebmt2)
#'     # transition matrix for competing risks
#'     tmat <- trans.comprisk(6,names=c("Relapse","GvHD","Bacterial","Viral","Fungal","Other"))
#'     # preparing long dataset
#'     ebmt2$stat1 <- as.numeric(ebmt2$status==1)
#'     ebmt2$stat2 <- as.numeric(ebmt2$status==2)
#'     ebmt2$stat3 <- as.numeric(ebmt2$status==3)
#'     ebmt2$stat4 <- as.numeric(ebmt2$status==4)
#'     ebmt2$stat5 <- as.numeric(ebmt2$status==5)
#'     ebmt2$stat6 <- as.numeric(ebmt2$status==6)
#'     covs <- c("dissub","match","tcd","year","age")
#'     ebmtlong <- msprep(time=c(NA,rep("time",6)),
#'                        stat=c(NA,paste("stat",1:6,sep="")),
#'                        data=ebmt2,keep=covs,trans=tmat)
#' 
#'     # The reduced rank 2 solution
#'     rr2 <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year+age,
#'                    data=ebmtlong, R=2)
#'     rr3$Alpha; rr3$Gamma; rr3$Beta; rr3$loglik
#'     # The reduced rank 3 solution
#'     rr3 <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year+age,
#'                    data=ebmtlong, R=3)
#'     rr3$Alpha; rr3$Gamma; rr3$Beta; rr3$loglik
#'     # The reduced rank 3 solution, with no reduction on age
#'     rr3 <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year, full=~age,
#'                    data=ebmtlong, R=3)
#'     rr3$Alpha; rr3$Gamma; rr3$Beta; rr3$loglik
#'     # The full rank solution
#'     fullrank <- redrank(Surv(Tstart,Tstop,status) ~ dissub+match+tcd+year+age,
#'                         data=ebmtlong, R=6)
#'     fullrank$Beta; fullrank$loglik
#' }
#' 
#' @export redrank
`redrank` <- function(redrank, full = ~1, data, R,
            strata = NULL, Gamma.start, method = "breslow", eps = 1e-5, print.level = 1)
{
    if (!inherits(data, "msdata"))
        stop("'data' must be an 'msdata' object")
    trans <- attr(data, "trans")
    # Now we need only the data
    data <- as.data.frame(data)

    # Get model matrices of reduced rank and full rank parts
    mmrr <- model.matrix(redrank, data=data)
    mmrr <- mmrr[,-1,drop=FALSE] # without intercept
    p <- ncol(mmrr)
    if (p==0) stop("Empty reduced rank part; please consider full model")
    mmrr <- data.frame(mmrr)
    covs <- names(mmrr)
    mmfull <- model.matrix(full, data=data)
    mmfull <- mmfull[,-1,drop=FALSE] # without intercept
    p2 <- ncol(mmfull)
    # Construct working data
    if (p2>0) {
        mmfull <- data.frame(mmfull)
        fullcovs <- names(mmfull)
        rrdata <- as.data.frame(data[,c("id","from","to","trans","Tstart","Tstop","time","status")])
        rrdata <- cbind(rrdata,mmrr,mmfull)
        cols <- 8 + (1:p)
        fullcols <- 8 + p + (1:p2)
    }
    else {
        rrdata <- as.data.frame(data[,c("id","from","to","trans","Tstart","Tstop","time","status")])
        rrdata <- cbind(rrdata,mmrr)
        cols <- 8 + (1:p)
        fullcols <- NULL        
    }

    # Get and store whether clock is forward or reset
    cx <- coxph(redrank, data=data)
    if (attr(cx$y, "type") == "counting") clock <- "forward" else if (attr(cx$y, "type") == "right") clock <- "reset" else stop("Surv object should be either of type 'counting' or 'right'")

    # Preparations for iterative algorithm
    trans2 <- to.trans2(trans)
    K <- nrow(trans2)
    if (!is.null(dimnames(trans)))
        tnames <- paste(trans2$fromname,"->",trans2$toname)
    else tnames <- as.character(1:K)
    ### add to the data set R replicates of columns with covariates Z_1...Z_p
    colsR <- matrix(0,R,p)
    for (r in 1:R) {
        ncd <- ncol(rrdata)
        rrdata <- cbind(rrdata,rrdata[,cols])
        colsR[r,] <- ((ncd+1):(ncd+p))
        names(rrdata)[((ncd+1):(ncd+p))] <- paste(covs, as.character(r), sep=".rr")
    }
    if (missing(Gamma.start)) Gamma.start <- matrix(rnorm(R*K),R,K)
    Gamma.iter <- Gamma.start
    iter <- 1
    prev.loglik <- 0
    loglik <- 100
    Delta <- loglik - prev.loglik

    while(abs(Delta) > eps) {
        if (print.level > 0) {
            cat("\n\nIteration", iter, "... \n")
            flush.console()
        }
        iter <- iter + 1
        attr(rrdata, "trans") <- trans
        class(rrdata) <- c("msdata", "data.frame")
        ms.it <- redrank.iter(data = rrdata, cols = cols, colsR = colsR,
            fullcols = fullcols, R = R, clock = clock, strata = strata, Gamma.iter = Gamma.iter,
            print.level = print.level - 1)
        Gamma.iter <- ms.it$Gamma
        if (print.level > 0) cat("\nGamma = ", round(Gamma.iter, 5))
        prev.loglik <- loglik
        loglik <- ms.it$loglik
        Delta <-  - loglik + prev.loglik
        if (print.level > 0) {
            cat("\nPrevious loglik =", round(prev.loglik, 5), ", present loglik =", round(loglik, 5),
                " Delta = ", round(Delta, 8))
            flush.console()
        }
    }
    Alpha <- ms.it$Alpha
    Gamma.final <- as.matrix(ms.it$Gamma)
    B <- Alpha %*% Gamma.final
    ### To make Alpha and Gamma unique, use singular value decomposition
    svd.B <- svd(B)
    Gamma.final <- (diag(sqrt(svd.B$d)) %*% t(svd.B$v))[1:R,]
    Gamma.final <- matrix(Gamma.final,R,K)
    rnames <- paste("r",1:R,sep="")
    dimnames(Gamma.final) <- list(rnames,tnames)
    Alpha <- (svd.B$u %*% diag(sqrt(svd.B$d)))[,1:R]
    if (R>1)
        norm.Alpha <- apply(Alpha^2,2, function(x) sqrt(sum(x)))
    else norm.Alpha <- sqrt(sum(Alpha^2))
    norm.Alpha.mat <- matrix(norm.Alpha,p,R,byrow=TRUE)
    Alpha <- Alpha/norm.Alpha.mat
    dimnames(Alpha) <- list(covs,rnames)
    AlphaX <-  as.matrix(rrdata[,cols]) %*% Alpha
    Gamma.final <- Gamma.final * matrix(norm.Alpha,R,K)
    dimnames(B) <- list(covs,tnames)
    return(list(Alpha = Alpha, Gamma = Gamma.final, Beta = B,
            Beta2 = ms.it$Beta2, cox.itr1 = ms.it$cox.itr1,
            AlphaX = AlphaX, niter = iter, df = R*(p+K-R), loglik = loglik))
}
