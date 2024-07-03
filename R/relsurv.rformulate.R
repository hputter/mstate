`rformulate.mstate` <- function (formula, data = parent.frame(), ratetable, na.action, 
                               rmap, int, centered, cause) 
{
  # Function rformulate taken from relsurv and adapted for use in mstate.
  
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m <- m[c(1, match(c("formula", "data", "cause"), names(m), 
                    nomatch = 0))]
  m[[1L]] <- quote(stats::model.frame)
  Terms <- if (missing(data)) 
    terms(formula, specials = c("strata", "ratetable"))
  else terms(formula, specials = c("strata", "ratetable"), 
             data = data)
  Term2 <- Terms
  rate <- attr(Terms, "specials")$ratetable
  if (length(rate) > 1) 
    stop("Can have only 1 ratetable() call in a formula")
  if (!missing(rmap)) {
    if (length(rate) > 0) 
      stop("cannot have both ratetable() in the formula and a rmap argument")
    rcall <- rmap
    if (!is.call(rcall) || rcall[[1]] != as.name("list")) 
      stop("Invalid rcall argument")
  }
  else if (length(rate) > 0) {
    stemp <- untangle.specials(Terms, "ratetable")
    rcall <- as.call(parse(text = stemp$var)[[1]])
    rcall[[1]] <- as.name("list")
    Term2 <- Term2[-stemp$terms]
  }
  else rcall <- NULL
  if (is.ratetable(ratetable)) {
    israte <- TRUE
    dimid <- names(dimnames(ratetable))
    if (is.null(dimid)) 
      dimid <- attr(ratetable, "dimid")
    else attr(ratetable, "dimid") <- dimid
    temp <- match(names(rcall)[-1], dimid)
    if (any(is.na(temp))) 
      stop("Variable not found in the ratetable:", (names(rcall))[is.na(temp)])
    if (any(!(dimid %in% names(rcall)))) {
      to.add <- dimid[!(dimid %in% names(rcall))]
      temp1 <- paste(text = paste(to.add, to.add, sep = "="), 
                     collapse = ",")
      if (is.null(rcall)) 
        rcall <- parse(text = paste("list(", temp1, 
                                    ")"))[[1]]
      else {
        temp2 <- deparse(rcall)
        rcall <- parse(text = paste("c(", temp2, ",list(", 
                                    temp1, "))"))[[1]]
      }
    }
  }
  else stop("invalid ratetable")
  newvar <- all.vars(rcall)
  if (length(newvar) > 0) {
    tform <- paste(paste(deparse(Term2), collapse = ""), 
                   paste(newvar, collapse = "+"), sep = "+")
    m$formula <- as.formula(tform, environment(Terms))
  }
  m <- eval(m, parent.frame())
  n <- nrow(m)
  if (n == 0) 
    stop("data set has 0 rows")
  Y <- stats::model.extract(m, "response")
  offset <- model.offset(m)
  if (length(offset) == 0) 
    offset <- rep(0, n)
  if (!is.Surv(Y)) 
    stop("Response must be a survival object")
  Y.surv <- Y
  if (attr(Y, "type") == "right") {
    type <- attr(Y, "type")
    status <- Y[, 2]
    Y <- Y[, 1]
    start <- rep(0, n)
    ncol0 <- 2
  }
  else if (attr(Y, "type") == "counting") {
    type <- attr(Y, "type")
    status <- Y[, 3]
    start <- Y[, 1]
    Y <- Y[, 2]
    ncol0 <- 3
  }
  else stop("Illegal response value")
  if (any(c(Y, start) < 0)) 
    stop("Negative follow up time")
  if (max(Y) < 30) 
    warning("The event times must be expressed in days! (Your max time in the data is less than 30 days) \n")
  rdata <- data.frame(eval(rcall, m), stringsAsFactors = TRUE)
  rtemp <- match.ratetable.mstate(rdata, ratetable)
  R <- rtemp$R
  cutpoints <- rtemp$cutpoints
  if (is.null(attr(ratetable, "factor"))) 
    attr(ratetable, "factor") <- (attr(ratetable, "type") == 
                                    1)
  attr(ratetable, "dimid") <- dimid
  rtorig <- attributes(ratetable)
  nrt <- length(rtorig$dimid)
  wh.age <- which(dimid == "age")
  wh.year <- which(dimid == "year")
  if (length(wh.age) > 0) {
    if (max(R[, wh.age]) < 150 & stats::median(diff(cutpoints[[wh.age]])) > 
        12) 
      warning("Age in the ratetable part of the formula must be expressed in days! \n (Your max age is less than 150 days) \n")
  }
  if (length(wh.year) > 0) {
    if (min(R[, wh.year]) > 1850 & max(R[, wh.year]) < 2020 & 
        inherits(cutpoints[[wh.year]], "rtdate"))
      warning("The calendar year must be one of the date classes (Date, date, POSIXt)\n (Your variable seems to be expressed in years) \n")
  }
  if (nrt != ncol(R)) {
    nonex <- which(is.na(match(rtorig$dimid, attributes(ratetable)$dimid)))
    for (it in nonex) {
      if (rtorig$type[it] != 1) 
        warning(paste("Variable ", rtorig$dimid[it], 
                      " is held fixed even though it changes in time in the population tables. \n (You may wish to set a value for each individual and not just one value for all)", 
                      sep = ""))
    }
  }
  strats <- attr(Term2, "specials")$strata
  if (length(strats)) {
    temp_str <- untangle.specials(Term2, "strata", 1)
    if (length(temp_str$vars) == 1) 
      strata.keep <- m[[temp_str$vars]]
    else strata.keep <- strata(m[, temp_str$vars], shortlabel = TRUE, 
                               sep = ",")
    Term2 <- Term2[-temp_str$terms]
  }
  else strata.keep <- factor(rep(1, n))
  if (!missing(cause)) 
    strata.keep <- factor(rep(1, n))
  attr(Term2, "intercept") <- 1
  X <- model.matrix(Term2, m)[, -1, drop = FALSE]
  mm <- ncol(X)
  if (mm > 0 && !missing(centered) && centered) {
    mvalue <- colMeans(X)
    X <- X - rep(mvalue, each = nrow(X))
  }
  else mvalue <- double(mm)
  cause <- stats::model.extract(m, "cause")
  if (is.null(cause)) 
    cause <- rep(2, nrow(m))
  keep <- Y > start
  if (!missing(int)) {
    int <- max(int)
    status[Y > int * 365.241] <- 0
    Y <- pmin(Y, int * 365.241)
    keep <- keep & (start < int * 365.241)
  }
  if (any(start > Y) | any(Y < 0)) 
    stop("Negative follow-up times")
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]
    Y <- Y[keep]
    start <- start[keep]
    status <- status[keep]
    R <- R[keep, , drop = FALSE]
    strata.keep <- strata.keep[keep]
    offset <- offset[keep]
    Y.surv <- Y.surv[keep, , drop = FALSE]
    cause <- cause[keep]
    n <- sum(keep)
    rdata <- rdata[keep, ]
  }
  temp <- R
  names(temp) <- paste0("X", 1:ncol(temp))
  data <- data.frame(start = start, Y = Y, stat = status, 
                     temp)
  if (mm != 0) 
    data <- cbind(data, X)
  attr(ratetable, "cutpoints") <- lapply(cutpoints, function(x) {
    if (inherits(x, "rtabledate")) 
      class(x) <- "date"
    x
  })
  out <- list(data = data, R = R, status = status, start = start, 
              Y = Y, X = as.data.frame(X), m = mm, n = n, type = type, 
              Y.surv = Y.surv, Terms = Terms, ratetable = ratetable, 
              offset = offset, formula = formula, cause = cause, mvalue = mvalue, 
              strata.keep = strata.keep)
  na.action <- attr(m, "na.action")
  if (length(na.action)) 
    out$na.action <- na.action
  out
}