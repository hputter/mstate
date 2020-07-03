require(lattice)

plot.MarkovTest <- function(x, y, what=c("states", "overall"), idx=NULL, quantiles=TRUE, qsup, states,
                            xlab, ylab, main, ...)
{
  # Not ready for publication, see definition of q95
  what <- match.arg(what)
  B <- dim(x$n_wb_trace)[1]
  ny <- length(y)
  if (missing(xlab)) xlab <- "Time"
  if (missing(ylab)) ylab <- "Test statistic"
  if (missing(main)) main <- ""
  if (what=="states") {
    # dfr <- x$zbar
    qualset <- x$qualset
    J <- length(qualset)
    dfr <- data.frame(time=rep(y, J), zbar=as.numeric(x$zbar), qualstate=rep(qualset, each=ny), ct=0)
    lwd <- 2
    lty <- 1
    col <- 1
    if (quantiles) {
      dfrl1 <- data.frame(time=rep(y, J), zbar=as.numeric(x$est_quant[1, , ]), qualstate=rep(qualset, each=ny), ct=1)
      dfru1 <- data.frame(time=rep(y, J), zbar=as.numeric(x$est_quant[2, , ]), qualstate=rep(qualset, each=ny), ct=3)
      dfr <- rbind(dfr, dfrl1, dfru1)
      lwd <- c(lwd, 2, 2)
      lty <- c(lty, 3, 3)
      col <- c(col, 1, 1)
    }
    if (!missing(qsup)) {
      if (qsup %in% 1:dim(x$b_stat_wb)[2]) {
        q95 <- apply(x$b_stat_wb[, qsup, ], 2, quantile, 0.95)
        print(q95)
        dfrl2 <- data.frame(time=rep(y, J), zbar=rep(-q95, each=ny), qualstate=rep(qualset, each=ny), ct=2)
        dfru2 <- data.frame(time=rep(y, J), zbar=rep(q95, each=ny), qualstate=rep(qualset, each=ny), ct=4)
        dfr <- rbind(dfr, dfrl2, dfru2)
        lwd <- c(lwd, 2, 2)
        lty <- c(lty, 3, 3)
        col <- c(col, 1, 1)
      }
    }
    if (!is.null(idx)) {
      idx <- intersect(1:B, idx)
      nB <- length(idx)
      if (nB > 0) {
        dfrb <- data.frame(time=rep(rep(y, J), each=nB),
                           zbar=as.numeric(x$n_wb_trace[idx, , ]),
                           qualstate=rep(qualset, each=ny*nB),
                           ct=rep(-idx, ny*J))
        dfr <- rbind(dfrb, dfr)
        lwd <- c(rep(0.5, nB), lwd)
        lty <- c(rep(1, nB), lty)
        col <- c(rep(8, nB), col)
      }
    }
    # print(dim(dfr))
    if (missing(states)) dfr$qualstate <- factor(dfr$qualstate)
    else dfr$qualstate <- factor(dfr$qualstate, levels=qualset, labels=states[qualset])
    xyplot(zbar ~ time | qualstate, data=dfr, groups=ct, lwd=lwd, type="l", col=col, lty=lty,
           xlab=xlab, ylab=ylab, main=main)
  }
  else if (what=="overall") {
    dfr <- data.frame(time=y, zbar=as.numeric(x$obs_chisq_trace), ct=0)
    lwd <- 2
    lty <- 1
    col <- 1
    if (quantiles) {
      
      dfru <- data.frame(time=y, zbar=apply(x$nch_wb_trace, 2, quantile, probs=0.95), ct=-1)
      dfr <- rbind(dfr, dfru)
      lwd <- c(2, lwd)
      lty <- c(3, lty)
      col <- c(1, col)
    }
    if (!is.null(idx)) {
      idx <- intersect(1:B, idx)
      nB <- length(idx)
      if (nB > 0) {
        dfrb <- data.frame(time=rep(y, each=nB),
                           zbar=as.numeric(x$nch_wb_trace[idx, ]),
                           ct=rep(idx, ny))
        dfr <- rbind(dfrb, dfr)
        lwd <- c(lwd, rep(0.5, nB))
        lty <- c(lty, rep(1, nB))
        col <- c(col, rep(8, nB))
      }
    }
    xyplot(zbar ~ time, data=dfr, groups=ct, lwd=lwd, type="l", col=col, lty=lty,
           xlab=xlab, ylab=ylab, main=main)
  }
}

