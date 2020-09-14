##***********************************##
## ggplot2 version of plot.probtrans ##
##         using geom_ribbons        ##
##***********************************##

# For data.table
.datatable.aware = TRUE

# Start function
ggplot.probtrans <- function(x,
                             from = 1,
                             type = "filled", 
                             ord,
                             cols,
                             xlab = "Time",
                             ylab = "Probability",
                             xlim,
                             ylim,
                             lwd, # linewidth
                             lty,
                             cex, # size of label
                             legend,
                             legend.pos = "right",
                             conf.int = 0.95, # 0 if no confidence intervals
                             conf.type = "log",
                             label) { # specify "annotate" for labels
                  
  # For R cmd check
  low <- upp <- state <- prob <- CI_low <- CI_upp <- NULL
               
  # Get names of states and state numbers
  if (missing(legend)) {
    state_names <- dimnames(x$trans)[[1]]
    
    # If supplied, check it is right length
  } else if (length(legend) != dim(x$trans)[1]) {
    stop(paste("'Legend' should be a vector of length ", dim(x$trans)[1]))
    
  } else state_names <- legend
  
  # Check stacking order
  if (missing(ord)) ord <- 1:length(state_names)
  if (missing(cex)) cex <- 8
  if (missing(lwd)) lwd <- 0.5
  
  # Check label is supplied for correct plot type
  if (!missing(label) & !(type %in% c("filled", "stacked")))
    stop("Labels only valid for filled/stacked plots!")
  
  # Create copy of pb object so original is NOT effected by references updating
  pb_copy <- data.table::copy(x)
  
  # Prepare dataframe 
  df_steps <- prep_probtrans_df(
    obj = pb_copy,
    from = from,
    ord = ord, 
    state_names = state_names,
    conf.int = conf.int,
    conf.type = conf.type
  )
  
  # Get bounds of plot
  if (missing(xlim)) xlim <- c(0, max(df_steps$time))
  if (missing(ylim)) ylim <- c(0, 1)

  # Arrange colours based on states plotted
  n_states_plotted <- length(unique(as.character(df_steps$state)))
  
  # check colours
  if (missing(cols)) cols <- viridis::viridis_pal()(n_states_plotted)
  if (length(cols) != n_states_plotted) 
    stop(paste0("Length of col should be ", n_states_plotted))
  
  # Default colours different if line plot
  if (type %in% c("separate", "single")) {
    
    # Extend colourbrewer Dark2 palette
    if (n_states_plotted <= 8) {
      full_pal <- RColorBrewer::brewer.pal(8, "Dark2")
      cols <- full_pal[1:n_states_plotted]
    } else {
      dark2_extended <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Dark2")
      )
      cols <- dark2_extended(n_states_plotted)
    }
  }
  
  # Check linetype
  if (missing(lty)) lty <- rep(1, n_states_plotted)
  if (length(lty) != n_states_plotted) 
    stop(paste0("Length of lty should be ", n_states_plotted))
  
  # Make different plot types
  if (type == "stacked") {
    
    p <- make_labelled_plot(
      df_steps = df_steps, 
      state_names = state_names,
      xlim = xlim,
      lwd = lwd,
      cex = cex
    ) + 
      
      # Stacked so we just remove fill colour
      ggplot2::scale_fill_manual(
        values = rep(NA, length(state_names))
      ) + 
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab) +
      ggplot2::coord_cartesian(expand = 0, xlim = xlim, ylim = ylim) +
      ggplot2::theme(legend.position = "none")
      
    
  } else if (type == "filled") {
    
    # add if for label_type
    if (missing(label)) {
      
      p <- ggplot2::ggplot(data = df_steps, 
                        ggplot2::aes(time, ymin = low, ymax = upp, 
                                     fill = state)) +
        ggplot2::geom_ribbon(col = "black", size = lwd,
                             na.rm = T) + # for rel surv
        ggplot2::guides(fill = ggplot2::guide_legend("State", reverse = T)) + 
        ggplot2::theme(legend.position = legend.pos) +
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0)  +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::scale_fill_manual(values = cols)
      
    } else {
      
      # Labelled plot
      p <- make_labelled_plot(
        df_steps = df_steps, 
        state_names = state_names,
        xlim = xlim,
        lwd = lwd,
        cex = cex
      ) + 
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::ylab(ylab) +
        ggplot2::xlab(xlab) +
        ggplot2::coord_cartesian(expand = 0, xlim = xlim, ylim = ylim) +
        ggplot2::theme(legend.position = "none")
      
    }
    
  } else if (type == "single") {
    
    # Colour of ribbon
    col_ribb <- "grey70"
    if (conf.type == "none") col_ribb <- NA
    
    p <- ggplot2::ggplot(data = df_steps,
                         ggplot2::aes(time, prob, 
                                      col = state, 
                                      linetype = state)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_low, 
                                        ymax = CI_upp,), 
                           alpha = 0.5, fill = col_ribb,
                           col = NA) +
      ggplot2::geom_line(size = lwd) + # colour boundaries
      ggplot2::guides(col = ggplot2::guide_legend("State", reverse = T)) +
      ggplot2::theme(legend.position = legend.pos) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::scale_colour_manual("State", values = cols) +
      ggplot2::scale_linetype_manual("State", values = lty) +
      ggplot2::guides(
        col = ggplot2::guide_legend("State", reverse = T),
        linetype = ggplot2::guide_legend("State", reverse = T)
      )
    
  } else if (type == "separate") {
    
    # Colour of ribbon
    col_ribb <- "grey70"
    if (conf.type == "none") col_ribb <- NA
    
    p <- ggplot2::ggplot(
      data = df_steps,
      ggplot2::aes(time, prob, 
                   ymin = CI_low, 
                   ymax = CI_upp,
                   col = state, 
                   linetype = state)
    ) +
      ggplot2::geom_ribbon(alpha = 0.5, fill = col_ribb,
                           col = NA, na.rm = T) +
      ggplot2::geom_line(size = lwd, na.rm = T) +
      ggplot2::facet_wrap(. ~ state) +
      ggplot2::guides(col = ggplot2::guide_legend("State", reverse = T)) +
      ggplot2::theme(legend.position = legend.pos) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::scale_colour_manual("State", values = cols) +
      ggplot2::scale_linetype_manual("State", values = lty) +
      ggplot2::guides(
        col = ggplot2::guide_legend("State", reverse = T),
        linetype = ggplot2::guide_legend("State", reverse = T)
      )
    
  } else stop("Pick a valid plot type!")
  
  return(p) 
}



# Helper function to prep data for plotting 
prep_probtrans_df <- function(obj,
                              from,
                              ord,
                              state_names,
                              conf.type = "log",
                              conf.int = 0.95) {
  
  # For R cmd check
  . <- state_num <- state <- se_state <- prob <- se <- cum_probs <-  NULL
  
  # Read in probtrans object
  df <- data.table::setDT(obj[[from]])
  
  # Subset states with at least one non zero probabily
  zero_prob_cols <- apply(df, 2, function(col) !all(col == 0))
  
  # Condition pstate cols
  pstate_cols <- grepl(x = names(df), pattern = "pstate|time")
  condition_pstate <- zero_prob_cols & pstate_cols  
  
  # Conditions se cols
  se_cols <- grepl(x = names(df), pattern = "se|time")
  condition_se <- zero_prob_cols & se_cols
  
  # Prepare 
  df_pstate <- df %>% 
    
    # Check for inf values, and exclude
    .[time != Inf, .SD, .SDcols = condition_pstate] %>%
    
    # Make long format
    data.table::melt.data.table(
      id.vars       = "time",
      variable.name = "state",
      value.name    = "prob"
    ) %>% 
    .[, state_num := as.numeric(
      gsub(x = state, pattern = "pstate", replacement = "")
    )]
  
  # Se's long
  df_se <- df %>% 
    
    # Check for infinite time values - exclude
    .[time != Inf, .SD, .SDcols = condition_se] %>%
    
    # Make long format
    data.table::melt.data.table(
      id.vars       = "time",
      variable.name = "se_state",
      value.name    = "se"
    ) %>% 
    
    .[, state_num := as.numeric(
      gsub(x = se_state, pattern = "se", replacement = "")
    )]
  
  
  # Put se's and pstates together
  df_long <- data.table::merge.data.table(
    x = df_pstate, 
    y = df_se, 
    by = c("time", "state_num")
  ) %>% 
    
    # Order and label factor with state names
    .[, state := factor(
      state_num, 
      levels = ord, 
      labels = state_names[ord]
    )] %>% 
    
    # Add CI for probabilities, do this on log
    .[, ':=' (
      CI_low = make_prob_confint(
        prob, se, conf.type, conf.int, bound = "low"
      ),
      CI_upp = make_prob_confint(
        prob, se, conf.type, conf.int, bound = "upp"
      )
    ), by = .(time, state)] %>% 
    
    # Sort by time and state, 
    .[order(time, state)] %>% 
    .[, cum_probs := cumsum(prob), by = time] %>% 
    
    # Compute upper and lower bounds of ribbons
    .[, ':=' (
      low = c(0, cum_probs[-length(cum_probs)]),
      upp = c(cum_probs[-length(cum_probs)], 1)
    ), by = time] 
  
  
  # Df with steps for ribbon
  df_steps <- rbind(df_long, 
                    df_long) %>%
    .[order(time)] %>% 
    
    # Shift time by 1, creating steps, equi to dplyr::lead(time, n = 1)
    .[, time := data.table::shift(
      time, 
      fill = NA, 
      n = 1, 
      type = "lead"
    ), by = state] %>% 
    
    .[!is.na(time)]
  
  return(df_steps) 
}


make_labelled_plot <- function(df_steps,
                               state_names,
                               xlim,
                               cex,
                               lwd) {
  
  # For R cmd check
  . <- state <- low <- upp <- mean_cprob <- label <- NULL
  
  # Max follow-up, depends on xlim
  time_vec <- df_steps[df_steps$time <= xlim[2], ]$time
  max_t <- time_vec[length(time_vec)]

  # Prep labels - at end of follow-up
  df_labels <- df_steps %>% 
    .[, ':=' (
      label = ifelse(time == max_t, as.character(state), NA_character_),
      mean_cprob = (low + upp) / 2
    )]  
  
  # Begin plots
  p <- df_labels %>% 
    
    # Remove last row to avoid dupl labels
    .[time <= max_t, .SD[-.N], by = state] %>% 
    ggplot2::ggplot(ggplot2::aes(time, mean_cprob)) + 
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = low, ymax = upp, fill = state),
      col = "black", size = lwd
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, x = time), 
      na.rm = T,
      hjust = 1, 
      size = cex
    ) 
    
  return(p)
}


make_prob_confint <- function(prob, 
                              se, 
                              conf.type = "log",
                              conf.int = 0.95, 
                              bound) {
  # Get critical value
  if (!is.null(conf.int)) {
    crit <- qnorm((1 - conf.int) / 2, lower.tail = FALSE)
  }

  if (conf.type == "log") {
    low <- exp(log(prob) - crit * se / prob)
    upp <-  exp(log(prob) + crit * se / prob)
      
  } else if (conf.type == "plain") {
    low <- prob - crit * se
    upp <- prob + crit * se
    
  } else {
    low <- prob
    upp <- prob
  }
  
  # Check no vals above 1 or below zero
  low <- ifelse(low < 0, 0, low)
  upp <- ifelse(upp > 1, 1, upp)
  
  # Return upper of lower bound
  if (bound == "upp") {
    return(upp) 
  } else return(low)
}
