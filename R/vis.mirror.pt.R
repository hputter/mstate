#' Mirror plot comparing two probtrans objects
#' 
#' A mirror plot for comparing two different \code{"probtrans"} objects. Useful
#' for comparing predicted probabilities for different levels of a covariate,
#' or for different subgroups at some prediction time horizon.
#'
#' @param x A list of two \code{"probtrans"} objects. 
#' The first element will be on the left of the mirror plot, 
#' and the second on the right
#' @param titles A character vector c("Title for left", "Title for right")
#' @param size_titles Numeric, size of the title text
#' @param breaks_x_left Numeric vector specifying axis breaks on the left plot
#' @param breaks_x_right Numeric vector specifying axis breaks on the right plot
#' @param xlab A title for the x-axis, default is "Time"
#' @param ylab A title for the y-axis, default is "Probability"
#' @param legend.pos Position of the legend, default is "right"
#' @param horizon Numeric, position denoting (in time) where to symmetrically 
#' mirror the plots. Default is maximum follow-up of from both plots.
#' 
#' @inheritParams plot.probtrans
#' 
#' @seealso \code{\link{plot.probtrans}}
#'
#' @return A ggplot2 object.
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
#' # 1. Prepare reference patient data - both CCR5 genotypes
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
#' # Mirror plot split at 10 years - see vignette for more details
#' vis.mirror.pt(
#' x = list(pt.WW, pt.WM),
#' titles = c("WW", "WM"),
#' horizon = 10
#' )
#' 
#' @export
vis.mirror.pt <- function(x, 
                          titles,
                          size_titles = 5,
                          horizon = NULL,
                          breaks_x_left,
                          breaks_x_right,
                          from = 1,
                          cols,
                          ord,
                          xlab = "Time",
                          ylab = "Probability",
                          legend.pos = "right") {
  
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.", call. = FALSE)
  }
  
  # Check the number of states
  tmat <- x[[1]][["trans"]]
  state_names <- dimnames(tmat)[[1]]
  n_states_plotted <- 1 + sum(!is.na(tmat[seq(from, nrow(tmat), by = 1), ]))
  
  if (missing(cols)) cols <- set_colours(n_states_plotted, type = "areas")
  if (missing(ord)) ord <- seq_along(state_names) 
  
  # First build plot and inherit colours from left
  p_left <- plot.probtrans(x[[1]], use.ggplot = TRUE, ord = ord, cols = cols, from = from)
  p_right <- plot.probtrans(x[[2]], use.ggplot = TRUE, ord = ord, cols = cols, from = from)
  build_p1 <- ggplot2::ggplot_build(p_left)

  # Inherit also xlims for symmetry - edit here for xlim.
  left_xlim <- p_left$coordinates$limits$x
  right_xlim <- p_right$coordinates$limits$x
  
  # Horizon for symmetry
  if (!is.null(horizon)) {
    
    if (horizon > max(left_xlim) | horizon > max(right_xlim)) 
      stop(paste0("Horizon should be smaller than both ", max(left_xlim), " and ", max(right_xlim)))
    left_xlim[which(left_xlim == max(left_xlim))] <- horizon
    right_xlim[which(right_xlim == max(right_xlim))] <- horizon
  }
  
  # In limits from plots
  dat_left <- p_left$data[time <= max(left_xlim)]
  dat_left[time == max(dat_left$time), ]$time <- max(left_xlim)
  
  dat_right <- p_right$data[time <= max(right_xlim)] 
  dat_right[time == max(dat_right$time), ]$time <- max(right_xlim)
  
  # Check if same ord was used for both plots
  if (any(levels(dat_left$state) != levels(dat_right$state)))
    stop("Argument 'ord' must be the same for both plots.")
  
  # Get maximum time and set x axis
  ylim <- c(0, 1)
  
  # Prep df, right side needs to lagged
  main <- prep_compare_df(
    dat_left = dat_left, 
    dat_right = dat_right
  )
  
  # Read-in objs
  dat_main <- main$df_compare
  max_t <- main$max_t
  diff <- main$diff
  
  # Prep labels
  breaks_size <- c(floor(max_t / 3), floor(max(dat_right$time) / 3))
  if (any(breaks_size == 0)) breaks_size <- c(round(max_t / 3, digits = 1), round(max(dat_right$time) / 3, digits = 1))
  if (missing(breaks_x_left)) breaks_x_left <- seq(from = 0, to = max_t, by = breaks_size[1])
  if (missing(breaks_x_right)) breaks_x_right <- seq(from = 0, to = max(dat_right$time), by = breaks_size[2])
  
  # If not missing check labels are within bounds
  if (any(breaks_x_left > max_t)) 
    stop(paste0("Max follow up on left side is ", round(max_t, 3), ", make sure all breaks are smaller!"))
  if (any(breaks_x_right > max(dat_right$time))) 
    stop(paste0("Max follow up on right side is ", round(max(dat_right$time), 3), ", make sure all breaks are smaller!"))
  breakos <- c(breaks_x_left, max(dat_main$time) - breaks_x_right)
  labos <- c(breaks_x_left, breaks_x_right)
  
  # Build basic plot
  p_main <-  ggplot2::ggplot(data = dat_main, ggplot2::aes(x = .data$time)) + 
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$low, ymax = .data$upp, fill = .data$state), 
      col = "black",
      na.rm = TRUE
    ) +
    # Add divider segment
    ggplot2::geom_segment(
      x = max_t, 
      xend = max_t,
      y = ylim[1], 
      yend = ylim[2], 
      size = 2
    ) +
    ggplot2::scale_x_continuous(xlab, breaks = breakos, labels = labos) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
  
  if (missing(titles)) {
    
    p <- p_main +
      ggplot2::theme(
        legend.position = legend.pos, 
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(), 
        axis.line = ggplot2::element_blank()
      ) +
      ggplot2::coord_cartesian(expand = 0, ylim = ylim) 
      
    
    return(p)
    
  } else {
    
    # Position of titles - divide by max time for grob
    pos_title_left <- (max_t / 2) 
    pos_title_right <- (max(dat_main$time) - max(dat_right$time) / 2)  

    titles_df <- data.frame(
      labs = c(titles[1], titles[2]),
      pos_x = c(pos_title_left, pos_title_right),
      pos_y = rep(1.05, 2)
    )
    
    base_p <- p_main + 
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(30.5, 5.5, 5.5, 5.5), "points"),
        legend.position = legend.pos,
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::coord_cartesian(clip = "off", expand = 0, ylim = ylim) +
      ggplot2::geom_text(
        data = titles_df, 
        ggplot2::aes(x = .data$pos_x, y = .data$pos_y, label = .data$labs),
        hjust = 0.5, 
        size = size_titles
      )
    
      
    return(base_p)
  }
}


prep_compare_df <- function(dat_left, dat_right) {
  
  # For data.table
  . <- time <- time_orig <- state <- side <- NULL
  
  # Get maximum times
  p1_maxt <- max(dat_left$time)
  p2_maxt <- max(dat_right$time)
  
  # Set max_t to left one in any case
  max_t <- p1_maxt
  
  # Check shift  
  diff <- ifelse(p2_maxt != p1_maxt, p2_maxt - p1_maxt, 0)
  
  # Get righthand side back to original df
  df_p2 <- data.table::copy(dat_right) 
  
  df_p2[, time := data.table::shift(
    x = time, 
    fill = NA,
    n = 1, 
    type = "lag"
  ), by = state] 
  
  df_p2 <- unique(df_p2[!is.na(time)])

  # Rescale time
  df_p2[, ':=' (
    time_orig = time,
    time = -time + 2 * max_t + (diff + 0.01)
  )]
    
  # Shift all
  cols <- c("time", "time_orig")  
  df_p2 <- rbind(df_p2, df_p2) 
  data.table::setorder(df_p2, time)
  df_p2[, (cols) := Map(
    f = data.table::shift,
    x = .SD, 
    fill = time[1L], 
    n = 1, 
    type = "lag"
  ), by = state, .SDcols = cols] 

  df_p2[, side := "right"]
  
  # Left sides
  df_p1 <- dat_left[, ':=' (
    time_orig = time,
    side = "left"
  )]
  
  # Bind both sides
  plot_df <- rbind(df_p1, df_p2[!is.na(time)])
  
  # Return a list of useful things
  res <- list("df_compare" = plot_df, "max_t" = max_t, "diff" = diff)
  return(res)
}

