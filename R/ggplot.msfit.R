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
                         lwd = 1) {
  
  # Extract transmat, and get transition labels
  trans_legend <- to.trans2(x$trans)
  
  # Extract df and append transition labels
  df <- x$Haz
  
  if (missing(legend)) {
    df$trans_name <- trans_legend$transname[match(df$trans, 
                                                  trans_legend$transno)]
  } else {
    
    if (length(legend) != nrow(trans_legend))
      stop(paste0("Legend is a char. vector of length ", nrow(trans_legend)))
    if (length(unique(legend)) < length(legend))
      stop("Pick unique labels for transitions!")
    
    # Match labels and transition numbers
    leg <- cbind.data.frame(transno = trans_legend$transno, transname = legend)
    df$trans_name <- leg$transname[match(df$trans, leg$transno)]
  }
  
  df$trans_name <- as.factor(df$trans_name)
  
  # Remove possible inf time values
  if (any(df$time == Inf)) df <- df[-which(df$time == Inf), ]
  
  # Checks
  n_states_plotted <- length(levels(df$trans_name))
  
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
  
  # Check lwd
  if (missing(lwd)) lwd <- 1
  
  # Check xlim and ylim
  if (missing(xlim)) xlim <- c(0, max(df$time))
  if (missing(ylim)) ylim <- c(0, max(df$Haz))
  
  # Start plots
  if (type == "single") {
    
    p <- ggplot2::ggplot(data = df, 
                         ggplot2::aes(.data$time, .data$Haz, col = .data$trans_name,
                             linetype = .data$trans_name)) +
      ggplot2::geom_step(size = lwd) +
      ggplot2::coord_cartesian(expand = 0, xlim = xlim, ylim = ylim) +
      ggplot2::ylab(ylab) +
      ggplot2::theme(legend.position = legend.pos) 
      
  } else if (type == "separate") {
    
    p <- ggplot2::ggplot(data = df, 
                         ggplot2::aes(.data$time, .data$Haz, 
                                      col = .data$trans_name,
                             linetype = .data$trans_name)) +
      ggplot2::geom_step(size = lwd) +
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

