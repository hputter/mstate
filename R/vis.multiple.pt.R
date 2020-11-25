#' Visualise multiple probtrans objects
#' 
#' Helper function allow to visualise state probabilities for 
#' different reference patients/covariates. Multiple \code{"probtrans"} objects
#' are thus needed.
#'
#' @param x A list of \code{"probtrans"} objects
#' @param from The starting state from which the probabilities are used to plot
#' Numeric, as in \code{plot.probtrans}
#' @param to (Numeric) destination state
#' @param labels Character vector labelling each element of x (e.g. label
#' for a reference patient) - so labels = c("Patient 1", "Patient 2")
#' @param legend.title Character - title of legend 
#' @inheritParams plot.probtrans
#'
#' @return A ggplot object.
#' 
#' @author Edouard F. Bonneville \email{e.f.bonneville@@lumc.nl}
#' 
#' @examples 
#' 
#' library(ggplot2)
#' 
#' data("aidssi")
#' head(aidssi)
#' si <- aidssi
#' 
#' # Prepare transition matrix
#' tmat <- trans.comprisk(2, names = c("event-free", "AIDS", "SI"))
#' 
#' # Run msprep
#' si$stat1 <- as.numeric(si$status == 1)
#' si$stat2 <- as.numeric(si$status == 2)
#' 
#' silong <- msprep(
#' time = c(NA, "time", "time"), 
#' status = c(NA, "stat1", "stat2"), 
#' data = si, keep = "ccr5", trans = tmat
#' )
#' 
#' # Run cox model
#' silong <- expand.covs(silong, "ccr5")
#' c1 <- coxph(Surv(time, status) ~ ccr5WM.1 + ccr5WM.2 + strata(trans),
#'             data = silong)
#'             
#' # 1. Prepare patient data - both CCR5 genotypes
#' WW <- data.frame(
#' ccr5WM.1 = c(0, 0), 
#' ccr5WM.2 = c(0, 0), 
#' trans = c(1, 2), 
#' strata = c(1, 2)
#' )
#' 
#' WM <- data.frame(
#' ccr5WM.1 = c(1, 0), 
#' ccr5WM.2 = c(0, 1),
#' trans = c(1, 2), 
#' strata = c(1, 2)
#' )
#' 
#' # 2. Make msfit objects
#' msf.WW <- msfit(c1, WW, trans = tmat)
#' msf.WM <- msfit(c1, WM, trans = tmat)
#' 
#' # 3. Make probtrans objects
#' pt.WW <- probtrans(msf.WW, predt = 0)
#' pt.WM <- probtrans(msf.WM, predt = 0)           
#' 
#' # Plot - see vignette for more details
#' vis.multiple.pt(
#' x = list(pt.WW, pt.WM), 
#' from = 1,
#' to = 2, 
#' conf.type = "log",
#' cols = c(1, 2),
#' labels = c("Pat WW", "Pat WM"),
#' legend.title = "Ref patients"
#' )
#' 
#' @export
vis.multiple.pt <- function(x, 
                            from = 1,
                            to,
                            xlab = "Time",
                            ylab,
                            xlim = NULL,
                            ylim = NULL,
                            cols,
                            lwd,
                            labels,
                            conf.int = 0.95,
                            conf.type = c("log", "plain", "none"),
                            legend.title) {
  
  # For cmd check
  . <- time <- state <- prob <- CI_low <- CI_upp <- PT <-  NULL
  
  # Check x are probtrans objects
  cond <- sapply(x, function(p) inherits(p, what = "probtrans"))
  if (!all(cond)) stop("x should be a list of probtrans objects")
  
  # Check conf
  conf.type <- match.arg(conf.type)
  
  # Check labels
  if (missing(labels)) labels <- paste0("pt_", 1:length(x))
  if (missing(legend.title)) legend.title <- "PT"
  if (missing(to)) stop("Please specify destination state in 'to'!")
  if (missing(lwd)) lwd <- 1
  
  # Check feasability of trans + get laveks
  tmat <- x[[1]]$trans
  if (is.na(tmat[from, to])) stop("This transition does not exist!")
  
  # Set to (using first pt)
  tmat_to <- dimnames(tmat)[["to"]]
  to <- tmat_to[to]
  
  # Set ylab
  if (missing(ylab)) ylab <- paste0(
    "Probability ", to, " from ", tmat_to[from]
  )
  
  # Number of pt
  n_pt <- length(x)
  
  # Check colours
  if (missing(cols)) {
    
    # Extend colourbrewer Dark2 palette
    if (n_pt <= 8) {
      full_pal <- RColorBrewer::brewer.pal(8, "Dark2")
      cols <- full_pal[1:n_pt]
    } else {
      dark2_extended <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Dark2")
      )
      cols <- dark2_extended(n_pt)
    }
    
  }
  
  # Prep data
  dfs_labelled <- lapply(1:length(x), function(i) {
    p <- plot(
      x[[i]], 
      type = "separate", 
      from = from, 
      conf.int = conf.int,
      conf.type = conf.type,
      use.ggplot = T
    )
    df <- p$data
    df[, PT := labels[i]]
    return(df)
  })
  
  # Plot
  p <- do.call(rbind, dfs_labelled) %>% 
    .[state == to] %>% 
    ggplot2::ggplot(ggplot2::aes(
      x = time, 
      y = prob, 
      ymin = CI_low, 
      ymax = CI_upp,
      group = PT
    )) +
    ggplot2::geom_ribbon(alpha = 0.5, fill = "grey70", col = NA) +
    ggplot2::geom_line(ggplot2::aes(col = PT), size = lwd) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::guides(colour = ggplot2::guide_legend(title = legend.title))
  
  return(p)
}