##***********************************##
## ggplot2 version of plot.Cuminc    ##
##         using geom_ribbons        ##
##***********************************##


prep_Cuminc_df <- function(x, 
                           shift_vars,
                           conf.type,
                           conf.int) {
  
  # For data.table warnings
  . <- prob <- se <- NULL
  
  # Get list of se cols and prob cols
  prob_cols <- grep(x = names(x), pattern = "^CI", value = T)
  se_cols <- grep(x = names(x), pattern = "^seCI", value = T)
  
  # Prepare long df
  df_long <- melt.data.table(
    data = data.table(x), 
    measure.vars = list(prob_cols, se_cols),
    #measure.vars = data.table:::patterns("^CI", "^seCI"), 
    value.name = c("prob", "se"), 
    variable.name = "state"
  ) %>% 
    .[, ':=' (
      CI_low = make_prob_confint(
        prob, se, conf.type, conf.int, bound = "low"
      ),
      CI_upp = make_prob_confint(
        prob, se, conf.type, conf.int, bound = "upp"
      )
    ), by = shift_vars] 
  
  # Shift for the ribbons
  df_steps <- rbind(df_long, 
                    df_long) %>%
    .[order(time)] %>% 
    
    # Shift time by 1, creating steps, equi to dplyr::lead(time, n = 1)
    .[, time := data.table::shift(
      time, 
      fill = NA, 
      n = 1, 
      type = "lead"
    ), by = shift_vars] %>% 
    
    .[!is.na(time)]
  
  return(df_steps)
} 


ggplot.Cuminc <- function(x,
                          xlab = "Time",
                          ylab = "Probability",
                          xlim,
                          ylim,
                          lty,
                          legend,
                          cols, 
                          conf.type = "log",
                          conf.int = 0.95,
                          legend.pos = "right",
                          facet = F) {
  
  # For data.table warnings
  state.grp <- state <- group <- NULL
  
  # Check whether group was specified
  if ("group" %in% names(x)) {
    shift_vars <- c("state", "group")
  } else shift_vars <- "state"
  
  # Check facet
  if (facet & length(shift_vars) == 1)
    stop("Cannot facet for Cuminc object without group specified")
  
  # Format data
  df_steps <- prep_Cuminc_df(
    x = x, 
    shift_vars = shift_vars, 
    conf.type = conf.type,
    conf.int = conf.int
  )
  
  # If group, add interaction variable for plotting 
  if (!facet & length(shift_vars) > 1) {
    df_steps[, state.grp := interaction(state, group)]
    grp <- rlang::sym("state.grp")
  } else grp <- rlang::sym("state")
  
  # Grps
  n_grps_plotted <- length(unique(df_steps[[grp]]))
  
  # Set up graphical parameters of plot
  if (missing(xlim)) xlim <- c(0, max(df_steps$time))
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lty)) lty <- rep(1, n_grps_plotted)
  if (missing(legend)) { legend <- levels(factor(df_steps[[grp]]))
    col_ribb <- NA
  } else col_ribb <- "grey70"
  
  # Colours
  if (missing(cols)) {
    
    # Extend colourbrewer Dark2 palette
    if (n_grps_plotted <= 8) {
      full_pal <- RColorBrewer::brewer.pal(8, "Dark2")
      cols <- full_pal[1:n_grps_plotted]
    } else {
      dark2_extended <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Dark2")
      )
      cols <- dark2_extended(n_grps_plotted)
    }
  }

  # Start plot
  p <- ggplot2::ggplot(
    data = df_steps,
    ggplot2::aes(x = .data$time,
                 y = .data$prob,
                 col = !!grp,
                 group = !!grp,
                 linetype = !!grp)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_upp),
      fill = col_ribb, 
      col = NA, 
      alpha = 0.5
    ) + 
    ggplot2::geom_line(size = 1) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0) +
    ggplot2::scale_linetype_manual(values = lty, labels = legend) +
    ggplot2::scale_color_manual(values = cols, labels = legend) +
    ggplot2::theme(legend.position = legend.pos) +
    ggplot2::labs(x = xlab, y = ylab)
    
  if (facet) {
    p <- p + ggplot2::facet_wrap(. ~ group)
  }
  
  return(p)
}

