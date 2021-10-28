#' Landmark Aalen-Johansen method
#' 
#' This function implements the landmark Aalen-Johansen method of Putter &
#' Spitoni (2016) for non-parametric estimation of transition probabilities in
#' non-Markov models.
#' 
#' 
#' @param msdata An \code{"msdata"} object, as for instance prepared by
#' \code{link{msprep}}
#' @param s The prediction time point s from which transition probabilities are
#' to be obtained
#' @param from Either a single state or a set of states in the state space
#' 1,...,S
#' @param method The method for calculating variances, as in
#' \code{\link{probtrans}}
#' @return A data frame containing estimates and associated standard errors of
#' the transition probabilities P(X(t)=k | X(s) in \code{from}) with \code{s}
#' and \code{from} the arguments of the function.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#' @references H. Putter and C. Spitoni (2016). Estimators of transition
#' probabilities in non-Markov multi-state models. Submitted.
#' @keywords survival
#' @examples
#' 
#' data(prothr)
#' tmat <- attr(prothr, "trans")
#' pr0 <- subset(prothr, treat=="Placebo")
#' attr(pr0, "trans") <- tmat
#' pr1 <- subset(prothr, treat=="Prednisone")
#' attr(pr1, "trans") <- tmat
#' c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=pr0)
#' c1 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=pr1)
#' msf0 <- msfit(c0, trans=tmat)
#' msf1 <- msfit(c1, trans=tmat)
#' # Comparison as in Figure 2 of Titman (2015)
#' # Aalen-Johansen
#' pt0 <- probtrans(msf0, predt=1000)[[2]]
#' pt1 <- probtrans(msf1, predt=1000)[[2]]
#' par(mfrow=c(1,2))
#' plot(pt0$time, pt0$pstate1, type="s", lwd=2, xlim=c(1000,4000), ylim=c(0,0.61),
#'      xlab="Time since randomisation (days)", ylab="Probability")
#' lines(pt1$time, pt1$pstate1, type="s", lwd=2, lty=3)
#' legend("topright", c("Placebo", "Prednisone"), lwd=2, lty=1:2, bty="n")
#' title(main="Aalen-Johansen")
#' # Landmark Aalen-Johansen
#' LMpt0 <- LMAJ(msdata=pr0, s=1000, from=2)
#' LMpt1 <- LMAJ(msdata=pr1, s=1000, from=2)
#' plot(LMpt0$time, LMpt0$pstate1, type="s", lwd=2, xlim=c(1000,4000), ylim=c(0,0.61),
#'      xlab="Time since randomisation (days)", ylab="Probability")
#' lines(LMpt1$time, LMpt1$pstate1, type="s", lwd=2, lty=3)
#' legend("topright", c("Placebo", "Prednisone"), lwd=2, lty=1:2, bty="n")
#' title(main="Landmark Aalen-Johansen")
#' 
#' @export LMAJ
LMAJ <- function(msdata, s, from, method=c("aalen", "greenwood"))
{
  tmat <- attr(msdata, "trans")
  if (is.null(tmat)) stop("msdata object should have a \"trans\" attribute")
  K <- nrow(tmat)
  if (any(is.na(match(from, 1:K)))) stop("from should be subset of 1:K with K number of states")
  xss <- xsect(msdata, s)
  infrom <- xss$id[xss$state %in% from]
  
  if (length(infrom) == 0) {
    msdata_from <- msdata[msdata$from == from, ]
    if (nrow(msdata_from) == 0) {
      stop(paste0("No transitions are made from state ", from, "!"))
    } else {
      first_entering <- round(min(msdata_from$Tstart), 4)
      stop_mssg <- paste0(
        "At landmark time s = ", s, 
        ", no individual has yet made it into state ", 
        from, ". The first individual enters at t = ", first_entering, "."
      )
      stop(stop_mssg)
    }
  } 
  
  msdatas <- cutLMms(msdata, LM=s)
  msdatasfrom <- msdatas[msdatas$id %in% infrom, ]
  msdatasfrom$trans <- factor(msdatasfrom$trans) # Feed to coxph factor (for easier matching)
  c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data=msdatasfrom)
  msf0 <- msfit(c0, trans=tmat)
  
  # Check if only subset of transitions left
  cond_subset_trans <- length(levels(msdatasfrom$trans)) < length(unique(msdata$trans))
  
  if (cond_subset_trans) {
    # Prepare original levels and re-coded ones
    recoded_levels <- levels(factor(as.numeric(msdatasfrom$trans)))
    orig_levels <- as.numeric(c0$xlevels$`strata(trans)`)

    # Match in the msf0 object accordingly
    msf0$Haz$trans <- orig_levels[match(as.character(msf0$Haz$trans), recoded_levels)]
    msf0$varHaz$trans1 <- orig_levels[match(as.character(msf0$varHaz$trans1), recoded_levels)]
    msf0$varHaz$trans2 <- orig_levels[match(as.character(msf0$varHaz$trans2), recoded_levels)]
  }
  
  # The warning is just for not being able to calculate variance at landmark time,
  # see probtrans.R, line 172
  pt0 <- probtrans(msf0, predt=s, method=method)[from]
  if (length(from) == 1)
    return(pt0[[1]])
  else {
    xsss <- xss[xss$state %in% from, ]
    xsss$state <- factor(xsss$state, levels=from)
    tbl <- table(xsss$state)
    p <- tbl / sum(tbl)
    varp <- (diag(p) - p %*% t(p)) / sum(tbl)
    # Relying on fact that all items from list have exactly the same size and structure
    # and that the sum of p equals 1; sorry for the double for-loop
    res <- tmp1 <- tmp2 <- 0
    for (j in 1:length(from)) {
      ptj <- pt0[[j]]
      res <- res + p[j] * ptj[, 1 + (1:K)]
      tmp2 <- tmp2 + p[j] * (1 - p[j]) * (ptj[, K+1 + (1:K)])^2
      for (k in 1:length(from)) {
        ptk <- pt0[[k]]
        tmp1 <- tmp1 + varp[j,k] * ptj[, 1 + (1:K)] * ptk[, 1 + (1:K)]
      }
    }
    ses <- sqrt(tmp1 + tmp2)
    # take ptj as template and insert estimates and SE's
    pt <- ptj
    pt[, 1 + (1:K)] <- res
    pt[, K+1 + (1:K)] <- ses
    return(pt)
  }
}
