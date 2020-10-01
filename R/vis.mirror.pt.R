#' Mirror plot comparing two probtrans objects
#' 
#' A mirror plot for comparing two different probtrans objects. Useful
#' for comparing predicted probabilities for different levels of a covariate,
#' or for different subgroups.
#'
#' @param x A list of two plots as returned by plot.probtrans with
#' use.ggplot = T. The first element will be on the left of the mirror plot,
#'  and the second on the right
#' @param titles A character vector c("Title for left", "Title for right")
#' @param size_titles Numeric, size of the title text
#' @param breaks_x_left Numeric vector specifying axis breaks on the left plot
#' @param breaks_x_right Numeric vector specifying axis breaks on the right plot
#' @param ylim Numeric vector, limits of the y-axis. Default is c(0, 1)
#' @param xlab A title for the x-axis, default is "Time"
#' @param ylab A title for the y-axis, default is "Probability"
#' @param legend.pos Position of the legend, default is "right"
#'
#' @return
#' @export
vis.mirror.pt <- function(x, 
                          titles,
                          size_titles = 15,
                          breaks_x_left,
                          breaks_x_right,
                          ylim,
                          xlab = "Time",
                          ylab = "Probability",
                          legend.pos = "right") {
  
   # First build plot and inherit colours from left
  p_left <- x[[1]]
  p_right <- x[[2]]
  build_p1 <- ggplot2::ggplot_build(p_left)
  fill_cols <- unique(build_p1$data[[1]]$fill)
  
  # Inherit also xlims for symmetry
  left_xlim <- p_left$coordinates$limits$x
  right_xlim <- p_right$coordinates$limits$x
  
  # In limits from plots
  dat_left <- p_left$data[time <= max(left_xlim)]
  dat_left[time == max(dat_left$time), ]$time <- max(left_xlim)
  
  dat_right <- p_right$data[time <= max(right_xlim)] 
  dat_right[time == max(dat_right$time), ]$time <- max(right_xlim)
  
  # Check if same ord was used for both plots
  if (any(levels(dat_left$state) != levels(dat_right$state)))
    stop("Argument 'ord' must be the same for both plots.")
  
  # Get maximum time and set x axis
  if (missing(ylim)) ylim <- c(0, 1)
  
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
  if (missing(breaks_x_left)) breaks_x_left <- seq(
    from = 0, to = max_t, by = floor(max_t / 3)
  )
  if (missing(breaks_x_right)) breaks_x_right <- seq(
    from = 0, to = max(dat_right$time), by = floor(max_t / 3)
  )
  
  # If not missing check labels are within bounds
  if (any(breaks_x_left > max_t))
    stop(paste0("Max follow up on left side is ", 
                round(max_t, 3),
                ", make sure all breaks are smaller!"))
  
  if (any(breaks_x_right > max(dat_right$time)))
    stop(paste0("Max follow up on right side is ", 
                round(max(dat_right$time), 3),
                ", make sure all breaks are smaller!"))
  
  breakos <- c(breaks_x_left, max(dat_main$time) - breaks_x_right)
  labos <- c(breaks_x_left, breaks_x_right)
  
  # Position of titles - divide by max time for grob
  pos_title_left <- (max_t / 2) / max(dat_main$time)
  pos_title_right <- (max(dat_main$time) - max(dat_right$time) / 2) / 
    max(dat_main$time)
  
  
  # Build basic plot
  p_main <-  dat_main %>% 
    ggplot2::ggplot(
      ggplot2::aes(x = .data$time, 
                   ymin = .data$low, 
                   ymax = .data$upp,
                   fill = .data$state)) + 
    ggplot2::geom_ribbon(col = "black", na.rm = T) +
    
    # Add divider segment
    ggplot2::geom_segment(
      x = max_t, 
      xend = max_t,
      y = ylim[1], 
      yend = ylim[2], 
      size = 2
    ) +
    ggplot2::coord_cartesian(expand = 0, ylim = ylim) +
    ggplot2::scale_x_continuous(
      xlab,
      breaks = breakos,
      labels = labos
    ) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_fill_manual(values = fill_cols) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = T))
  
  if (missing(titles)) {
    
    p <- p_main +
      ggplot2::theme(
        legend.position = legend.pos, 
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(), 
        axis.line = ggplot2::element_blank()
      )
    
    return(p)
    
  } else {
    
    base_p <- p_main + 
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(30.5, 5.5, 5.5, 5.5), "points"),
        legend.position = legend.pos,
        panel.grid = ggplot2::element_blank()
      )
    
    # See https://cran.r-project.org/web/packages/gridExtra/vignettes/gtable.html
    # https://stackoverflow.com/questions/21997715/add-ggplot-annotation-outside-the-panel-or-two-titles
    
    # Add titles using gtable
    g <- ggplot2::ggplotGrob(base_p) %>% 
      gtable::gtable_add_grob(
        grid::grobTree(
          grid::textGrob(
            label = titles[1], 
            x = pos_title_left, 
            hjust = 0.5,
            gp = grid::gpar(fontsize = size_titles)
          ), 
          grid::textGrob(
            label = titles[2], 
            x = pos_title_right, 
            hjust = 0.5,
            gp = grid::gpar(fontsize = size_titles)
          )
        ), t = 1, l = 5
      )
    
    grid::grid.draw(g)
    return(invisible(g))
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
  if (p2_maxt != p1_maxt) {
    diff <- p2_maxt - p1_maxt
  } else diff <- 0
  
  # Get righthand side back to original df
  df_p2 <- data.table::copy(dat_right) %>% 
    .[, time := data.table::shift(
      time, 
      fill = NA,
      n = 1, 
      type = "lag"
    ), by = state] %>% 
    .[!is.na(time)] %>% 
    unique() %>% 
    
    # Rescale time
    .[, ':=' (
      time_orig = time,
      time = -time + 2 * max_t + (diff + 0.01)
    )]
  
  cols <- c("time", "time_orig")  
  
  df_p2 <- rbind(df_p2, df_p2) %>% 
    .[order(time)] %>% 
    .[, (cols) := Map(
      data.table::shift,
      .SD, 
      fill = time[1L], 
      n = 1, 
      type = "lag"
    ), by = state, .SDcols = cols] %>% 
    .[!is.na(time)]  %>% 
    .[, side := "right"]
  
  df_p1 <- dat_left[, ':=' (
    time_orig = time,
    side = "left"
  )]
  
  plot_df <- rbind(df_p1, df_p2)
  
  # Return a list of useful things
  res <- list("df_compare" = plot_df,
              "max_t" = max_t,
              "diff" = diff)
}

