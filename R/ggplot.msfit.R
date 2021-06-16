##*******************************##
## ggplot2 version of plot.msfit ##
##         using geom_ribbons    ##
##*******************************##


ggplot.msfit <- function(x,
                         type = "single", # or "separate"
                         xlab = "Time", 
                         ylab = "Cumulative hazard", 
                         xlim,
                         ylim,
                         scale_type = "fixed", # see scales arg of facet_wrap
                         legend,
                         legend.pos = "right",
                         lty,
                         cols,
                         lwd = 1,
                         conf.int = 0.95, 
                         conf.type = "log") {
  
  # Extract transmat, and get transition labels
  trans_legend <- to.trans2(x$trans)
  
  # Extract df and append transition labels
  df <- x$Haz
  
  if (missing(legend)) {
    df$trans_name <- trans_legend$transname[match(df$trans, trans_legend$transno)]
  } else {
    
    if (length(legend) != nrow(trans_legend)) stop(paste0("Legend is a char. vector of length ", nrow(trans_legend)))
    if (length(unique(legend)) < length(legend)) stop("Pick unique labels for transitions!")
    
    # Match labels and transition numbers
    leg <- cbind.data.frame(transno = trans_legend$transno, transname = legend)
    df$trans_name <- leg$transname[match(df$trans, leg$transno)]
  }
  
  df$trans_name <- factor(df$trans_name, levels = unique(df$trans_name))
  
  # Remove possible inf time values
  if (any(df$time == Inf)) df <- df[-which(df$time == Inf), ]
  n_states_plotted <- length(levels(df$trans_name))
  
  # Merge with variance part
  df_steps <- prep_msfit_df(
    df_haz = df, 
    df_var = x$varHaz, 
    conf.int = conf.int, 
    conf.type = conf.type
  )
  
  # Check colours
  if (missing(cols)) {
    cols <- set_colours(n_states_plotted, type = "lines")
  } else if (length(cols) != n_states_plotted) {
    stop(paste0("Length of col should be ", n_states_plotted))
  }
  
  # Check linetype
  if (missing(lty)) {
    lty <- rep(1, n_states_plotted)
  } else if (length(lty) != n_states_plotted) {
    stop(paste0("Length of lty should be ", n_states_plotted))
  }
  
  # Check lwd, xlim and ylim
  if (missing(lwd)) lwd <- 1
  if (missing(xlim)) xlim <- c(0, max(df_steps$time, na.rm = TRUE))
  if (missing(ylim)) ylim <- c(0, max(df_steps$CI_upp, na.rm = TRUE))
  col_ribb <- ifelse(conf.type == "none", NA, "grey70")
  
  # Start plots
  if (type == "single") {
    
    p <- ggplot2::ggplot(
      data = df_steps,
      ggplot2::aes(
        x = .data$time, 
        y = .data$Haz, 
        col = .data$trans_name, 
        linetype = .data$trans_name,
        group = .data$trans_name
      )
    ) + 
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_upp, group = .data$trans_name), 
        alpha = 0.5, 
        fill = col_ribb,
        col = NA,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(size = lwd) +
      ggplot2::coord_cartesian(expand = 0, xlim = xlim, ylim = ylim) +
      ggplot2::ylab(ylab) +
      ggplot2::theme(legend.position = legend.pos) 
      
  } else if (type == "separate") {
    
    p <- ggplot2::ggplot(
      data = df_steps, 
      ggplot2::aes(
        x = .data$time, 
        y = .data$Haz, 
        col = .data$trans_name,
        linetype = .data$trans_name
      )
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_upp), 
        alpha = 0.5, 
        fill = col_ribb,
        col = NA,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(size = lwd) +
      ggplot2::ylab(ylab) +
      ggplot2::coord_cartesian(expand = 0, xlim = xlim) +
      ggplot2::facet_wrap(. ~ trans_name, scales = scale_type) +
      ggplot2::theme(legend.position = "none") 
     
  } else stop("Invalid plot type!")
 
  p <- p +
    ggplot2::scale_colour_manual("Transition", values = cols) +
    ggplot2::scale_linetype_manual("Transition", values = lty) +
    ggplot2::xlab(xlab) 
  
  return(p) 
}


make_haz_confint <- function(haz, 
                             se, 
                             conf.type = "log",
                             conf.int = 0.95, 
                             bound) {
  # Get critical value
  crit <- if (!is.null(conf.int)) {
    qnorm((1 - conf.int) / 2, lower.tail = FALSE)
  } else 0
  
  if (conf.type == "log") {
    low <- exp(log(haz) - crit * se / haz)
    upp <-  exp(log(haz) + crit * se / haz)
    
  } else if (conf.type == "plain") {
    low <- haz - crit * se
    upp <- haz + crit * se
    
  } else {
    low <- haz
    upp <- haz
  }
  
  # Return upper of lower bound
  if (bound == "upp") {
    return(upp) 
  } else return(low)
}


prep_msfit_df <- function(df_haz, df_var, conf.type, conf.int) {
  
  varHaz <- Haz <- se <- trans_name <- NULL
  
  # Merge df with variances to the one with hazards
  df_varHaz <- df_var[df_var$trans1 == df_var$trans2, ]
  df_varHaz$trans <- df_varHaz$trans1 
  df_long <- data.table::data.table(merge(x = df_haz, y = df_varHaz, by = c("time", "trans")))
  
  # Make confidence intervals
  df_long[, "se" := sqrt(varHaz)]
  df_long[, ':=' (
    CI_low = make_haz_confint(Haz, se, conf.type, conf.int, bound = "low"),
    CI_upp = make_haz_confint(Haz, se, conf.type, conf.int, bound = "upp")
  )]
  
  # Df with steps for ribbon
  df_steps <- rbind(df_long, df_long) 
  data.table::setorder(df_steps, time, trans_name)
  
  # Shift time by 1, creating steps, equi to dplyr::lead(time, n = 1)
  df_steps[, time := data.table::shift(
    x = time, 
    fill = NA, 
    n = 1, 
    type = "lead"
  ), by = trans_name] 
  
  df_steps <- df_steps[!is.na(time)]
  
  return(df_steps)
}
