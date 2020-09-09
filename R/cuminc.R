#' Calculate nonparametric cumulative incidence functions and associated
#' standard errors
#' 
#' This function computes nonparametric cumulative incidence functions and
#' associated standard errors for each value of a group variable.
#' 
#' The estimated cumulative incidences are as described in Putter, Fiocco &
#' Geskus (2007); the standard errors are the square roots of the Greenwood
#' variance estimators, see eg. Andersen, Borgan, Gill & Keiding (1993), de
#' Wreede, Fiocco & Putter (2009), and they correspond to the variances in eg.
#' Marubini & Valsecchi (1995). In case of no censoring, the estimated
#' cumulative incidences and variances reduce to simple binomial frequencies
#' and their variances.
#' 
#' @aliases Cuminc print.Cuminc plot.Cuminc summary.Cuminc
#' @param time Either 1) a numeric vector containing the failure times or 2) a
#' string containing the column name indicating these failure times
#' @param status Either 1) a numeric, factor or character vector containing the
#' failure codes or 2) a string containing the column name indicating these
#' failure codes
#' @param data When appropriate, a data frame containing \code{time},
#' \code{status} and/or \code{group} variables
#' @param group Optionally, name of column in data indicating a grouping
#' variable; cumulative incidence functions are calculated for each value or
#' level of \code{group}. If missing no groups are considered
#' @param failcodes A vector indicating which values of \code{status} are
#' considered as different causes of failure; other values of \code{status} are
#' considered as censorings. If missing and \code{status} is numeric, it is
#' assumed that 0 is censoring and all other values indicate failcodes; if
#' missing and \code{status} is character or factor, then it is assumed that
#' each of the levels/values of \code{status} is a cause of failure
#' @param na.status One of \code{"remove"} (default) or \code{"extra"},
#' indicating whether subjects with missing cause of failure should be removed
#' or whether missing cause of failure should be treated as a separate cause of
#' failure
#' @param variance Logical value, indicating whether the standard errors of the
#' cumulative incidences should be output (\code{TRUE}, the default) or not
#' @param x Object of class \code{"Cuminc"} to be printed or plotted
#' @param object Object of class \code{"Cuminc"} to be summarized
#' @param \dots Further arguments to plot or print method
#' @return An object of class \code{"Cuminc"}, which is a data frame containing
#' the estimated failure-free probabilities and cumulative incidences and their
#' standard errors. The names of the dataframe are \code{time}, \code{Surv},
#' \code{seSurv}, and \code{cuminc} and \code{secuminc} followed by the values
#' or levels of the \code{failcodes}. If \code{group} was specified, a
#' \code{group} variable is included with the same name and values/levels as
#' the original grouping variable, and with estimated cumulative incidences
#' (SE) for each value/level of \code{group}.
#' 
#' Cuminc is now simply a wrapper around survfit of the survival package with
#' type=\code{"mstate"}, only maintained for backward compatibility. The
#' survfit object is kept as attribute (\code{attr("survfit")}), and the print,
#' plot and summary functions are simply print, plot and summary applied to the
#' survfit object. Subsetting the \code{"Cuminc"} object results in subsetting
#' the data frame, not in subsetting the survfit object.
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @references Andersen PK, Borgan O, Gill RD, Keiding N (1993).
#' \emph{Statistical Models Based on Counting Processes}. Springer, New York.
#' 
#' Marubini E, Valsecchi MG (1995). \emph{Analysing Survival Data from Clinical
#' Trials and Observational Studies}. Wiley, New York.
#' 
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: Competing
#' risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
#' 2389--2430.
#' 
#' de Wreede L, Fiocco M, Putter H (2009). The mstate package for estimation
#' and prediction in non- and semi-parametric multi-state models. Submitted.
#' \url{http://www.msbi.nl/multistate}.
#' @keywords survival
#' @examples
#' 
#' ### These data were used in Putter, Fiocco & Geskus (2007)
#' data(aidssi)
#' ci <- Cuminc(time=aidssi$time, status=aidssi$status)
#' head(ci); tail(ci)
#' ci <- Cuminc(time="time", status="status", data=aidssi, group="ccr5")
#' head(ci); tail(ci)
#' 
#' ### Some fake data
#' fake <- data.frame(surv=c(seq(2,10,by=2),seq(1,13,by=3),seq(1,9,by=2),seq(1,13,by=3)),
#'                     stat=rep(0:3,5),Tstage=c(1:4,rep(1:4,rep(4,4))))
#' fake$stat[fake$stat==0 & fake$Tstage==2] <- 3
#' fake$stat[fake$stat==3 & fake$Tstage==1] <- 2
#' fake
#' Cuminc(time="surv", status="stat", data=fake)
#' # If we remove all entries with status=0,
#' # we should get binomial sample probabilities and corresponding SEs
#' fake0 <- fake[fake$stat!=0,]
#' Cuminc(time="surv", status="stat", data=fake0)
#' 
#' @export Cuminc
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
  if (length(status) != length(time))
    stop("lengths of time and status do not match")
  n <- length(time)
  
  if (missing(group)) {
    tmp <- data.frame(time=time, status=status)
    # Just call survfit with status argument a factor
    tmp$statuscr <- factor(tmp$status)
    sf <- survfit(Surv(time, statuscr) ~ 1, data=tmp)
    tt <- sf$time
    CIs <- sf$pstate
    ses <- sf$std.err
    res <- cbind(tt, CIs, ses)
    res <- as.data.frame(res)
    names(res) <- c("time", "Surv", paste("CI", sf$states[-1], sep="."),
                    "seSurv", paste("seCI", sf$states[-1], sep="."))
    res <- res[!duplicated(res$Surv), ]
    rownames(res) <- 1:nrow(res)
    class(res) <- c("Cuminc", "data.frame")
    attr(res, "survfit") <- sf
  } else {
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
    group <- as.factor(group)
    if (is.null(groupname)) groupname <- "group"
    tmp <- data.frame(time=time, status=status, group=group)
    # Call survfit with status argument a factor
    tmp$statuscr <- factor(tmp$status)
    sf <- survfit(Surv(time, statuscr) ~ group, data=tmp)
    tt <- sf$time
    CIs <- sf$pstate
    ses <- sf$std.err
    res <- cbind(tt, CIs, ses)
    res <- as.data.frame(res)
    names(res) <- c("time", "Surv", paste("CI", sf$states[-1], sep="."),
                    "seSurv", paste("seCI", sf$states[-1], sep="."))
    group <- rep(levels(as.factor(tmp$group)), sf$strata)
    res <- cbind(data.frame(group=group), res)
    res <- res[!duplicated(res$Surv), ]
    rownames(res) <- 1:nrow(res)
    class(res) <- c("Cuminc", "data.frame")
    attr(res, "survfit") <- sf
  }
  return(res)
}
