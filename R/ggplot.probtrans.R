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
  } else if (length(legend) != dim(x$trans)[1]) {
    stop(paste("'Legend' should be a vector of length ", dim(x$trans)[1]))
  } else state_names <- legend
  
  # Check stacking order
  if (missing(ord)) ord <- seq_along(state_names) 
  if (missing(cex)) cex <- 8
  if (missing(lwd)) lwd <- 0.5
  if (!missing(label) & !(type %in% c("filled", "stacked"))) stop("Labels only valid for filled/stacked plots!")
  
  # Create copy of pb object so original is NOT affected by references updating
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
  
  # Set graphical parameters
  if (missing(xlim)) xlim <- c(0, max(df_steps$time, na.rm = TRUE))
  if (missing(ylim)) ylim <- c(0, 1)
  n_states_plotted <- length(unique(as.character(df_steps$state)))
  if (missing(cols)) cols <- set_colours(n_states_plotted, type = "areas")
  if (length(cols) != n_states_plotted) stop(paste0("Length of col should be ", n_states_plotted))
  if (type %in% c("separate", "single")) cols <- set_colours(n_states_plotted, type = "lines")
  if (missing(lty)) lty <- rep(1, n_states_plotted)
  if (length(lty) != n_states_plotted) stop(paste0("Length of lty should be ", n_states_plotted))
  
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
      ggplot2::scale_fill_manual(values = rep("transparent", length(state_names))) + 
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab) +
      ggplot2::coord_cartesian(expand = 0, xlim = xlim, ylim = ylim) +
      ggplot2::theme(legend.position = "none")
      
    
  } else if (type == "filled") {
    
    # add if for label_type
    if (missing(label)) {
      
      p <- ggplot2::ggplot(
        data = df_steps, 
        ggplot2::aes(
          x = .data$time, 
          ymin = .data$low, 
          ymax = .data$upp, 
          fill = .data$state
        )
      ) +
        ggplot2::geom_ribbon(col = "black", size = lwd, na.rm = TRUE) + # for rel surv
        ggplot2::guides(fill = ggplot2::guide_legend("State", reverse = TRUE)) + 
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
    col_ribb <- ifelse(conf.type == "none", NA, "grey70")
    
    p <- ggplot2::ggplot(
      data = df_steps,
      ggplot2::aes(
        x = .data$time, 
        y = .data$prob, 
        col = .data$state, 
        linetype = .data$state
      )
    ) 
    
    if (sum(grepl(pattern = "CI_low|CI_upp", x = names(df_steps))) > 0) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_upp), 
        alpha = 0.5, 
        fill = col_ribb,
        col = NA,
        na.rm = TRUE
      )
    }
    
    p <- p +
      ggplot2::geom_line(size = lwd) + # colour boundaries
      ggplot2::guides(col = ggplot2::guide_legend("State", reverse = TRUE)) +
      ggplot2::theme(legend.position = legend.pos) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::scale_colour_manual("State", values = cols) +
      ggplot2::scale_linetype_manual("State", values = lty) +
      ggplot2::guides(
        col = ggplot2::guide_legend("State", reverse = TRUE),
        linetype = ggplot2::guide_legend("State", reverse = TRUE)
      )
    
  } else if (type == "separate") {
    
    # Colour of ribbon
    col_ribb <- ifelse(conf.type == "none", NA, "grey70")
    
    p <- ggplot2::ggplot(
      data = df_steps,
      ggplot2::aes(
        x = .data$time,
        y = .data$prob, 
        col = .data$state, 
        linetype = .data$state
      )
    ) 
    
    if (sum(grepl(pattern = "CI_low|CI_upp", x = names(df_steps)) > 0)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_upp), 
        alpha = 0.5, 
        fill = col_ribb,
        col = NA,
        na.rm = TRUE
      )
    }
    
    p <- p +
      ggplot2::geom_line(size = lwd, na.rm = TRUE) +
      ggplot2::facet_wrap(. ~ state) +
      ggplot2::guides(col = ggplot2::guide_legend("State", reverse = TRUE)) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = 0) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::scale_colour_manual("State", values = cols) +
      ggplot2::scale_linetype_manual("State", values = lty) +
      ggplot2::guides(
        col = ggplot2::guide_legend("State", reverse = TRUE),
        linetype = ggplot2::guide_legend("State", reverse = TRUE)
      ) +
      # Facet titles are already the labels
      ggplot2::theme(legend.position = "none") 
      
    
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
  df <- data.table::data.table(obj[[from]])
  
  # Subset states with at least one non zero probability
  zero_prob_cols <- apply(df, 2, function(col) !all(col == 0))
  
  # Condition pstate cols
  pstate_cols <- grepl(x = names(df), pattern = "pstate|time")
  condition_pstate <- zero_prob_cols & pstate_cols  

  # Prepare 
  df_pstate <- data.table::melt.data.table(
    data = df[time != Inf, .SD, .SDcols = condition_pstate],
    id.vars = "time",
    variable.name = "state",
    value.name = "prob"
  ) 
  df_pstate[, state_num := as.numeric(gsub(x = state, pattern = "pstate", replacement = ""))]
  df_long <- df_pstate

  # Check if standard errors were computed (only time)
  se_cols <- grepl(x = names(df), pattern = "se|time")
  if (sum(se_cols) > 1) {
    condition_se <- zero_prob_cols & se_cols
    
    # Se's long
    df_se <- data.table::melt.data.table(
      data = df[time != Inf, .SD, .SDcols = condition_se],
      id.vars = "time",
      variable.name = "se_state",
      value.name = "se"
    ) 
    df_se[, state_num := as.numeric(gsub(x = se_state, pattern = "se", replacement = ""))]
    
    # Put se's and pstates together
    df_long <- data.table::merge.data.table(
      x = df_pstate, 
      y = df_se, 
      by = c("time", "state_num")
    )
    
    # Add CI for probabilities, do this on log
    df_long[, ':=' (
      CI_low = make_prob_confint(prob, se, conf.type, conf.int, bound = "low"),
      CI_upp = make_prob_confint(prob, se, conf.type, conf.int, bound = "upp")
    )]  #, by = .(time, state)] %>% 
  }
  
  # Order and label factor with state names
  df_long[, state := factor(
    x = state_num, 
    levels = ord, 
    labels = state_names[ord]
  )]
  
  data.table::setorder(df_long, time, state)
  df_long[, cum_probs := cumsum(prob), by = time]
  
  # Compute upper and lower bounds of ribbons
  df_long[, ':=' (
    low = c(0, cum_probs[-length(cum_probs)]),
    upp = c(cum_probs[-length(cum_probs)], 1)
  ), by = time] 
  
  # Df with steps for ribbon
  df_steps <- rbind(df_long, df_long) 
  data.table::setorder(df_steps, time)
    
  # Shift time by 1, creating steps, equi to dplyr::lead(time, n = 1)
  df_steps[, time := data.table::shift(
    x = time, 
    fill = NA, 
    n = 1, 
    type = "lead"
  ), by = state] 
 
  return(df_steps[!is.na(time)]) 
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
  df_labels <- df_steps[, ':=' (
    label = ifelse(time == max_t, as.character(state), NA_character_),
    mean_cprob = (low + upp) / 2
  )]  
  
  # Begin plots
  p <- ggplot2::ggplot(
    data = df_labels[time <= max_t, .SD[-.N], by = state],
    ggplot2::aes(x = .data$time, y = .data$mean_cprob)
  ) + 
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$low, ymax = .data$upp, fill = .data$state),
      col = "black", size = lwd
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label, x = .data$time), 
      na.rm = T,
      hjust = 1, 
      size = cex
    ) 
    
  return(p)
}


make_prob_confint <- function(prob, 
                              se, 
                              conf.type = c("log", "plain", "none"),
                              conf.int = 0.95, 
                              bound = c("low", "upp")) {
  
  conf.type <- match.arg(conf.type)
  bound <- match.arg(bound)
  
  # Get critical value
  crit <- ifelse(!is.null(conf.int), qnorm((1 - conf.int) / 2, lower.tail = FALSE), 0L)
  direc <- ifelse(bound == "low", -1L, 1L)
  
  bound <- switch(
    conf.type,
    log = exp(log(prob) + direc * crit * se / prob),
    plain = prob + direc * crit * se,
    none = prob
  )

  # Bound limits in [0, 1], and no CIs for estimates with no SEs
  bound[bound < 0] <- 0
  bound[bound > 1] <- 1
  bound[se == 0] <- prob[se == 0]

  return(bound)
}
