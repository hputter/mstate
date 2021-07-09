`match.ratetable.mstate` <- function (R, ratetable) 
{
  # Function match.ratetable taken from survival and adapted for use in mstate.
  
  datecheck <- function(x) inherits(x, c("Date", "POSIXt", 
                                         "date", "chron", "rtabledate"))
  if (!is.ratetable(ratetable)) 
    stop("Invalid rate table")
  dimid <- names(dimnames(ratetable))
  if (is.null(dimid)) 
    dimid <- attr(ratetable, "dimid")
  datecut <- sapply(attr(ratetable, "cutpoints"), datecheck)
  rtype <- attr(ratetable, "type")
  if (is.null(rtype)) {
    temp <- attr(ratetable, "factor")
    rtype <- 1 * (temp == 1) + ifelse(datecut, 3, 2) * (temp == 
                                                          0) + 4 * (temp > 1)
  }
  if (is.matrix(R)) {
    attR <- attributes(R)
    attributes(R) <- attR["dim"]
    Rnames <- attR$dimnames[[2]]
    isDate <- attR[["isDate"]]
    levlist <- attR[["levlist"]]
  }
  else {
    Rnames <- names(R)
    levlist <- lapply(R, levels)
    isDate <- sapply(R, datecheck)
  }
  ord <- match(dimid, Rnames)
  if (any(is.na(ord))) 
    stop(paste("Argument '", dimid[is.na(ord)], "' needed by the ratetable was not found in the data", 
               sep = ""))
  if (any(duplicated(ord))) 
    stop("A ratetable argument appears twice in the data")
  R <- R[, ord, drop = FALSE]
  levlist <- levlist[ord]
  isDate <- isDate[ord]
  dtemp <- dimnames(ratetable)
  if (any((rtype < 3) & isDate)) {
    indx <- which(rtype < 3 & isDate)
    stop(paste("Data has a date type variable, but the reference", 
               "ratetable is not a date variable:", paste(dimid[indx], 
                                                          collapse = " ")))
  }
  if (any((rtype > 2) & !isDate)) {
    indx <- which(rtype > 2 & !isDate)
  }
  for (i in (1:ncol(R))) {
    if (rtype[i] > 2) 
      R[, i] <- ratetableDate(R[, i])
    if (length(levlist[[i]]) > 0) {
      if (rtype[i] != 1) 
        stop(paste("for this ratetable,", dimid[i], 
                   "must be a continuous variable"))
      temp <- charmatch(casefold(levlist[[i]]), casefold(dtemp[[i]]))
      if (any(is.na(temp))) 
        stop(paste("Levels do not match for ratetable() variable", 
                   dimid[i]))
      if (any(temp == 0)) 
        stop(paste("Non-unique ratetable match for variable", 
                   dimid[i]))
      R[, i] <- temp[as.numeric(R[, i])]
    }
    else {
      R[, i] <- unclass(R[, i])
      if (rtype[i] == 1) {
        temp <- R[, i]
        if (any(floor(temp) != temp) || any(temp <= 0) || 
            max(temp) > length(dtemp[[i]])) 
          stop(paste("The variable", dimid[i], 
                     "is out of range"))
      }
    }
  }
  R <- as.matrix(R)
  summ <- function(R){
    x <- c(format(round(min(R[, 1])/365.241, 1)), 
           format(round(max(R[,1])/365.241, 1)), sum(R[, 3] == 1), sum(R[, 3] == 2))
    x2 <- as.character(as.Date(c(min(R[, 2]), max(R[, 2])), origin=as.Date('1970-01-01')))
    paste("  age ranges from", x[1], "to", x[2], "years\n", " male:", 
          x[3], " female:", x[4], "\n", " date of entry from", 
          x2[1], "to", x2[2], "\n")
  }
  cutpoints <- lapply(attr(ratetable, "cutpoints"), ratetableDate)
  if (is.null(summ)) 
    list(R = R, cutpoints = cutpoints)
  else list(R = R, cutpoints = cutpoints, summ = summ(R))
}