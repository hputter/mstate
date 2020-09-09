#' Log-rank based test for the validity of the Markov assumption
#' 
#' Log-rank based test for the validity of the Markov assumption
#' 
#' Function MarkovTest performs the log-rank test described in Titman & Putter
#' (2020). Function optimal_weights_matrix implements the optimal weighting for
#' the state-specific trace. Function optimal_weights_multiple implements the
#' optimal weighting for the chi-squared trace.
#' 
#' @aliases MarkovTest optimal_weights_multiple optimal_weights_matrix
#' @param data Multi-state data in \code{msdata} format. Should also contain
#' (dummy codings of) the relevant covariates; no factors allowed
#' @param id Column name in \code{data} containing subject id
#' @param formula Right-hand side of the formula. If NULL will fit with no
#' covariates (formula="1" will also work), offset terms can also be specified.
#' @param transition Transition number of the transition to be tested (in the
#' transition matrix as attribute to \code{data})
#' @param grid Grid of time points at which to compute the statistic
#' @param B Number of wild bootstrap replications to perform
#' @param fn A list of summary functions to be applied to the individual zbar
#' traces (or a list of lists)
#' @param fn2 A list of summary functions to be applied to the overall
#' chi-squared trace
#' @param dist Distribution of wild bootstrap random weights, either "poisson"
#' for centred Poisson (default), or "normal" for standard normal
#' @param min_time The minimum time for calculating optimal weights
#' @param other_weights Other (than optimal) weights can be specified here
#' @return MarkovTest returns an object of class "MarkovTest", which is a list
#' with the following items: \item{orig_stat}{Summary statistic for each of the
#' starting states} \item{orig_ch_stat}{Overall chi-squared summary statistic}
#' \item{p_stat_wb}{P-values corresponding to each of the summary statistics
#' for each starting state} \item{p_ch_stat_wb}{P-values for overall
#' chi-squared summary statistic} \item{b_stat_wb}{Bootstrap summary statistics
#' for each of the starting states} \item{zbar}{Individual traces for each of
#' the starting states} \item{nobs_grid}{The number of events after time s for
#' each s in the grid} \item{Nsub}{Number of patients who are ever at risk of
#' the transition of interest} \item{est_quant}{Pointwise 2.5 and 97.5 quantile
#' limits for each of the traces} \item{obs_chisq_trace}{Trace of the
#' chi-squared statistic} \item{nch_wb_trace}{Individual values of the
#' chi-squared statistic trace for the wild bootstrap samples}
#' \item{n_wb_trace}{Individual values of the log-rank z statistic traces for
#' the wild bootstrap samples} \item{est_cov}{Estimated covariance matrix
#' between the log-rank statistics at each grid point} \item{transition}{The
#' transition number tested} \item{from}{The from state of the transition
#' tested} \item{to}{The to state of the transition tested} \item{B}{The number
#' of wild bootstrap replications} \item{dist}{The distribution used in the
#' wild bootstrap} \item{qualset}{Set of qualifying states corresponding to the
#' components of the above traces} \item{coxfit}{Fitted coxph object}
#' \item{fn}{List of functions applied to state-specific trace} \item{fn2}{List
#' of functions applied to overall trace}
#' @author Andrew Titman \email{a.titman@@lancaster.ac.uk}, transported to
#' mstate by Hein Putter \email{H.Putter@@lumc.nl}
#' @references Titman AC, Putter H (2020). General tests of the Markov property
#' in multi-state models. \emph{Biostatistics} To appear.
#' @keywords univar
#' @examples
#' 
#' # Example provided by the prothrombin data
#' data("prothr")
#' # Apply Markov test to grid of monthly time points over the first 7.5 years
#' year <- 365.25
#' month <- year / 12
#' grid <- month * (1 : 90)
#' # Markov test for transition 1 (wild bootstrap based on 25 replications, 1000 recommended)
#' MT <- MarkovTest(prothr, id = "id", transition = 1,
#'                  grid = grid, B = 25)
#' 
#' # Plot traces
#' plot(MT, grid, what="states", idx=1:10, states=rownames(attr(prothr, "trans")),
#'      xlab="Days since randomisation", ylab="Log-rank test statistic",
#'      main="Transition Normal -> Low")
#' plot(MT, grid,what="overall", idx=1:10,
#'      xlab="Days since randomisation", ylab="Chi-square test statistic",
#'      main="Transition Normal -> Low")
#' 
#' # Example using optimal weights and adjustment for covariates
#' oweights_fun <-
#'   optimal_weights_matrix(prothr, id = "id", grid=grid, transition = 1,
#'                          other_weights=list(
#'                            function(x) mean(abs(x),na.rm=TRUE),
#'                            function(x) max(abs(x),na.rm=TRUE)))
#' 
#' oweights_chi <- optimal_weights_multiple(prothr, id = "id", grid=grid, transition = 1)
#' 
#' # Formula in MarkovTest only works for continuous covariates and dummy coded variables
#' # No factors allowed
#' prothr$prednisone <- as.numeric(prothr$treat == "Prednisone")
#' MT <- MarkovTest(prothr, id = "id", 
#'                  formula = "prednisone",
#'                  transition = 1,
#'                  grid = grid, B = 25,
#'                  fn = oweights_fun,
#'                  fn2 = list(
#'                    function(x) weighted.mean(x, w=oweights_chi, na.rm=TRUE),
#'                    function(x) mean(x, na.rm=TRUE),
#'                    function(x) max(x, na.rm=TRUE)))
#' 
MarkovTest <- function(data, id, formula = NULL, transition, grid, 
                       B = 1000,
                       fn = list(function(x) mean(abs(x), na.rm = TRUE)),
                       fn2 = list(function(x) mean(x, na.rm = TRUE)),
                       dist = c("poisson", "normal")) {
  
  dist <- match.arg(dist)
  if (missing(id)) id <- "id"
  # Remove "empty" lines in the data
  wh <- which(data$Tstop <= data$Tstart)
  if (length(wh)>0)
  {
    warning(length(wh), " lines with Tstart <= Tstop, have been removed before applying tests!")
    data <- data[-wh, ]
  }

  # Convert data to etm data
  # Make sure to retain all covariates (possibly way to many) in msdata (needed in formula perhaps)
  mtch <- match(c("id", "from", "to", "trans", "Tstart", "Tstop", "status"), names(data))
  covcols <- 1:ncol(data)
  covcols <- covcols[!covcols %in% mtch]
  ncovs <- length(covcols)
  
  trans <- attr(data, "trans")
  etmdata <- msdata2etm(data, id)
  if (ncovs > 0) etmdata <- msdata2etm(data, id, names(data)[covcols])
  trans2 <- to.trans2(trans)
  tfrom <- trans2$from[trans2$transno == transition]
  tto <- trans2$to[trans2$transno == transition]
  
  # Determine qualifying set
  qualset <- c(tfrom, which(trans2Q(trans)[, tfrom] > 0))
  qualset <- sort(unique(qualset))  # for circular models, tfrom is included twice

  # Functions
  if (!is.list(fn)) 
    fn <- list(fn)  # coerce to be list if a single function is provided
  if (is.list(fn) & is.function(fn[[1]])) {
    # coerce to be a list of lists, by repeating the same list each time
    tempfn <- list()
    for (i in 1:length(qualset)) tempfn[[i]] <- fn
    fn <- tempfn
  }
  if (!is.list(fn2)) 
    fn2 <- list(fn2)  # coerce to be list if a single function is provided
  
  # Establish the relevant patients who ever enter tfrom
  relpat <- sort(unique(etmdata$id[etmdata$from == tfrom]))
  
  rdata <- etmdata[etmdata$from == tfrom, ]  # only need time periods in the relevant state...
  rdata$status <- 1 * (rdata$to == tto)
  if (!is.null(formula)) {
    form <- as.formula(paste("Surv(entry, exit, status) ~ ", formula, 
                             sep = ""))
    progfit <- coxph(form, data = rdata)
    if (length(progfit$coefficients) > 0) {
      Zmat <- as.matrix(rdata[, match(names(progfit$coefficients), 
                                      names(rdata))])
      Ncov <- dim(Zmat)[2]
    } else {
      Ncov <- 0
    }
    if (!is.null(progfit$offset)) {
      offset <- progfit$offset
    } else {
      offset <- rep(0, dim(rdata)[1])
    }
  } else {
    Ncov <- 0
    offset <- rep(0, dim(rdata)[1])
    progfit <- NULL
  }
  
  # Minimal data, change names
  progdat <- rdata[, match(c("id", "entry", "exit", "status"), names(rdata))]
  names(progdat) <- c("id", "T0", "T1", "D")
  
  nobs_grid <- sapply(grid, function(x) sum(progdat$D[progdat$T1 > x]))
  
  # Have the extra dimension of indexes
  index_gM <- array(0, c(length(relpat), length(grid), length(qualset)))
  for (indx in 1:length(qualset)) {
    qualstate <- qualset[indx]
    index_g <- sapply(grid, function(y) sapply(relpat, function(x)
      which(etmdata$entry < y & etmdata$exit >= y & etmdata$id == x)))
    index_g <- array(1 * (etmdata$from[sapply(index_g, function(y)
      ifelse(length(y) > 0, y, 
             dim(etmdata)[1] + 1))] == qualstate), c(length(relpat), length(grid)))
    index_g[is.na(index_g)] <- 0
    index_gM[, , indx] <- index_g
  }

  # Need a separate Z3mat for each group as well...
  Z3mat <- index_gM[match(progdat$id, relpat), , , drop = FALSE]
  
  N1 <- dim(progdat)[1]
  
  if (Ncov > 0) {
    LP <- c(Zmat %*% progfit$coefficients) + offset
  } else {
    LP <- rep(0, N1) + offset
  }
  S0 <- sapply(1:N1, function(x) sum(exp(LP) * (progdat$T0 < progdat$T1[x] & 
                                                  progdat$T1 >= progdat$T1[x])))
  
  incr <- progdat$D / S0
  cumhaz <- approxfun(c(0, sort(unique(progdat$T1)), Inf),
                      c(0, cumsum(tapply(incr, progdat$T1, sum)), sum(incr)),
                      method = "constant")
  resid_mat <- sapply(grid, function(x) progdat$D * (progdat$T1 > x) - exp(LP) *
                        (cumhaz(pmax(x, progdat$T1)) - cumhaz(pmax(x, progdat$T0))))
  
  # Have a separate trace for each qualifying state...
  obs_trace <- array(0, c(length(grid), length(qualset)))
  for (indx in 1:(length(qualset))) {
    obs_trace[, indx] <- sapply(1:length(grid), function(k)
      sum(resid_mat[, k] * Z3mat[, k, indx] * (progdat$T1 > grid[k])))
  }
  
  nqstate <- length(qualset)
  if (Ncov > 0) 
    Ifish <- progfit$var
  
  N1 <- dim(progdat)[1]
  if (Ncov > 0) 
    Zbar0 <- array(0, c(N1, Ncov))
  
  Zbar <- array(0, c(N1, length(grid), nqstate))
  for (i in 1:N1) {
    x <- i
    if (Ncov > 0) {
      for (j in 1:Ncov) {
        Zbar0[i, j] <- sum(Zmat[, j] * exp(LP) *
                             (progdat$T0 < progdat$T1[x] & progdat$T1 >= progdat$T1[x])) /
          sum(exp(LP) * (progdat$T0 < progdat$T1[x] & progdat$T1 >= progdat$T1[x]))
      }
    }
    for (j in 1:length(grid)) {
      for (k in 1:nqstate)
        Zbar[i, j, k] <- sum(Z3mat[, j, k] * exp(LP) * 
                               (progdat$T0 < progdat$T1[x] & progdat$T1 >= progdat$T1[x])) /
          sum(exp(LP) * (progdat$T0 < progdat$T1[x] & progdat$T1 >= progdat$T1[x]))
    }
  }
  
  NAe <- incr
  
  if (Ncov > 0) {
    Hmat <- array(0, c(length(grid), Ncov, nqstate))
    for (j in 1:Ncov) {
      # for (k in 1:nqstate) Hmat[,j,k] <- sapply(1:length(grid),function(y)
      # sum(sapply(1:N1,function(x) sum(exp(LP) *Zmat[,j]* (Z3mat[x,y,k] -
      # Zbar[x,y,k]) * NAe[x] * (progdat$T0[x] > grid[y] & progdat$T1[x] <=
      # progdat$T1)))))
      for (k in 1:nqstate) Hmat[, j, k] <- sapply(1:length(grid), function(y) 
        sum(sapply(1:N1, function(x)
          sum(exp(LP[x]) * ((Zmat[x, j] - Zbar0[, j]) *
                              (Z3mat[x, y, k] - Zbar[, y, k])) * NAe *
                (progdat$T1[x] > grid[y]) * (progdat$T1 > progdat$T0[x] & progdat$T1 <= progdat$T1[x])))))
    }
  }
  
  if (Ncov > 0) {
    multiplier <- array(0, dim(Hmat))
    for (k in 1:nqstate) multiplier[, , k] <- Hmat[, , k] %*% Ifish
    est_cov <- array(0, c(length(grid), nqstate, nqstate))
    for (indx1 in 1:nqstate) {
      for (indx2 in (indx1):nqstate) {
        est_var <- sapply(1:length(grid), function(k)
          sum(sapply(1:N1, function(v) 
            sum(((Z3mat[v, k, indx1] - Zbar[, k, indx1]) * 
                   (progdat$T1 > grid[k]) - c(multiplier[k, , indx1, drop = FALSE] %*%
                                                t(Zmat[v, ] - Zbar0))) *
                  ((Z3mat[v, k, indx2] - Zbar[, k, indx1]) * (progdat$T1 > grid[k]) - 
                     c(multiplier[k, , indx2, drop = FALSE] %*% t(Zmat[v, ] - Zbar0))) *
                  exp(LP[v]) * (progdat$T0[v] < progdat$T1 & progdat$T1[v] >= progdat$T1) * NAe))))
        est_cov[, indx1, indx2] <- est_cov[, indx2, indx1] <- est_var
      }
    }
    
  } else {
    est_cov <- array(0, c(length(grid), nqstate, nqstate))
    for (indx1 in 1:nqstate) {
      for (indx2 in (indx1):nqstate) {
        est_var <- sapply(1:length(grid), function(k)
          sum(sapply(1:N1, function(v)
            sum((Z3mat[v, k, indx1] - Zbar[, k, indx1]) * (Z3mat[v, k, indx2] - Zbar[, k, indx2]) *
                  exp(LP[v]) * (progdat$T1 > grid[k] & progdat$T0[v] < progdat$T1 & progdat$T1[v] >= progdat$T1) * NAe))))
        est_cov[, indx1, indx2] <- est_cov[, indx2, indx1] <- est_var
      }
    }
  }
  
  # First obtain the individually normalized traces...
  est_var <- obs_norm_trace <- array(0, c(length(grid), nqstate))
  for (k in 1:nqstate) {
    est_var[, k] <- est_cov[cbind(1:length(grid), k, k)]
    # This should be the same as before...
    obs_norm_trace[, k] <- obs_trace[, k] / sqrt(est_var[, k] + 1 * (est_var[, k] == 0))
  }
  # Find singular matrices
  obs_chisq_trace <- rep(0, length(grid))
  for (k in 1:length(grid)) {
    sol <- tryCatch(solve(est_cov[k, -1, -1]), error = function(e)
      return(diag(0, nqstate - 1)))
    obs_chisq_trace[k] <- (obs_trace[k, -1]) %*% sol %*%
      (obs_trace[k, -1]) # do something about singular matrices...
  }
  
  ############## 
  
  n_wb_trace <- wb_trace0 <- wb_trace <- array(0, c(B, length(grid), nqstate))
  nch_wb_trace <- array(0, c(B, length(grid)))
  for (wb in 1:B) {
    if (dist == "poisson") {
      G <- rpois(dim(progdat)[1], 1) - 1
    } else if (dist == "normal") {
      G <- rnorm(dim(progdat)[1], 0, 1)
    } else stop("argument dist should be poisson or normal")
    trace0 <- array(0, c(length(grid), nqstate))
    for (k in 1:nqstate) {
      trace0[, k] <- apply(sapply(1:length(grid), function(x)
        progdat$D * (Z3mat[, x, k] - Zbar[, x, k]) * (progdat$T1 > grid[x]) * G), 2, sum)
      if (Ncov > 0) {
        Imul <- sapply(1:Ncov, function(x)
          sum(progdat$D * (Zmat[, x] - Zbar0[, x]) * G))
        trace1 <- (Hmat[, , k] %*% Ifish %*% Imul)[, 1]
      } else {
        trace1 <- 0
      }
      wb_trace[wb, , k] <- trace0[, k] - trace1
      n_wb_trace[wb, , k] <- wb_trace[wb, , k]/sqrt(est_var[, k] + 
                                                      1 * (est_var[, k] == 0))
      for (w in 1:length(grid)) {
        sol <- tryCatch(solve(est_cov[w, -1, -1]), error = function(e)
          return(diag(0, nqstate - 1)))
        nch_wb_trace[wb, w] <- (wb_trace[wb, w, -1]) %*% sol %*% 
          (wb_trace[wb, w, -1]) # do something about singular matrices...
      }
    }
  }
  
  # Need to have one of these per nqstate
  NS <- length(fn[[1]])
  
  orig_stat <- array(sapply(1:nqstate, function(y)
    sapply(fn[[y]], function(g) g(obs_norm_trace[, y]))), c(NS, nqstate))
  orig_ch_stat <- sapply(fn2, function(g) g(obs_chisq_trace))
  
  p_stat_wb <- array(0, c(NS, nqstate))
  wb_stat <- array(0, c(B, NS, nqstate))
  for (k in 1:nqstate) {
    wb_stat[, , k] <- array(t(apply(n_wb_trace[, , k, drop = FALSE], 
                                    1, function(x)
                                      sapply(fn[[k]], function(g) g(x)))), c(B, NS))
    p_stat_wb[, k] <- sapply(1:NS, function(x) mean(wb_stat[, x, k] > orig_stat[x, k]))
  }
  est_quant <- array(0, c(2, length(grid), nqstate))
  for (k in 1:nqstate)
    est_quant[, , k] <- apply(n_wb_trace[, , k, drop = FALSE], 2,
                              quantile, c(0.025, 0.975), na.rm = TRUE)
  NS2 <- length(fn2)
  p_ch_stat_wb <- rep(0, NS2)
  wb_ch_stat <- array(t(apply(nch_wb_trace, 1, function(x)
    sapply(fn2, function(g) g(x)))), c(B, NS2))
  p_ch_stat_wb <- sapply(1:NS2, function(x) mean(wb_ch_stat[, x] > orig_ch_stat[x]))
  # Is a question whether should use Nsub as number of subjects or number
  # of spells within the state
  
  MTres <- list(orig_stat = orig_stat, orig_ch_stat = orig_ch_stat, p_stat_wb = p_stat_wb, 
                p_ch_stat_wb = p_ch_stat_wb, b_stat_wb = wb_stat, zbar = obs_norm_trace, 
                nobs_grid = nobs_grid, Nsub = length(relpat), est_quant = est_quant, 
                obs_chisq_trace = obs_chisq_trace, nch_wb_trace = nch_wb_trace, 
                n_wb_trace = n_wb_trace, est_cov = est_cov, transition = transition,
                from = tfrom, to = tto, B = B, dist = dist,
                qualset = qualset, coxfit = progfit, fn = fn, fn2 = fn2)
  
  class(MTres) <- c("MarkovTest")
  return(MTres)
}

optimal_weights_multiple <- function(data, id, grid, transition, min_time = 0)
{
  # Convert data to etm data
  trans <- attr(data, "trans")
  etmdata <- msdata2etm(data, id)
  trans2 <- to.trans2(trans)
  from <- trans2$from[trans2$transno == transition]
  to <- trans2$to[trans2$transno == transition]

  numbers <- sapply(grid, function(x)
    table(factor(etmdata$from)[(etmdata$entry <= x & etmdata$exit > x)]))
  subevent <- sapply(grid, function(x)
    sum(etmdata$from == from & etmdata$to == to & etmdata$exit > x))
  tnumbers <- apply(numbers, 2, sum)
  weights <- sapply(1:dim(numbers)[1], function(x)
    subevent * numbers[x, ] * (tnumbers - numbers[x, ])/tnumbers^2)
  weights[is.nan(weights)] <- 0
  weight <- apply(weights, 1, max)
  weight * diff(c(min_time, grid))
}

optimal_weights_matrix <- function(data, id, grid, transition, min_time = 0, 
                                   other_weights = NULL)
{
  # Convert data to etm data
  trans <- attr(data, "trans")
  etmdata <- msdata2etm(data, id)
  trans2 <- to.trans2(trans)
  from <- trans2$from[trans2$transno == transition]
  to <- trans2$to[trans2$transno == transition]

  numbers <- sapply(grid, function(x)
    table(factor(etmdata$from)[(etmdata$entry <= x & etmdata$exit > x)]))
  subevent <- sapply(grid, function(x)
    sum(etmdata$from == from & etmdata$to == to & etmdata$exit > x))
  tnumbers <- apply(numbers, 2, sum)
  weights <- sapply(1:dim(numbers)[1], function(x)
    sqrt(subevent * numbers[x, ] * (tnumbers - numbers[x, ]))/tnumbers)
  weights[is.nan(weights)] <- 0
  fn_list <- list()
  for (i in 1:dim(numbers)[1]) {
    # Take into account the distance between grids
    val <- weights[, i] * diff(c(min_time, grid))
    fn_list[[i]] <- list(fn = function(x)
      weighted.mean(abs(x), w = val, na.rm = TRUE))
    if (!is.null(other_weights)) {
      nother <- length(other_weights)
      fn_list[[i]][2:(nother + 1)] <- other_weights
    }
  }
  # Store the weights as an attribute
  attr(fn_list, "weights") <- weights
  fn_list
}
