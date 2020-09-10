#' Visualise multiple probtrans objects
#' 
#' Description here.
#'
#' @param x A list of plots as returned by plot(pt, use_ggplot = T)
#' @param to Character; destination state
#' @param labels Character vector labelling each element of x (e.g. label
#' for a reference patient) - so labels = c("Patient 1", "Patient 2")
#'
#' @return A plot visualising 
#' @export
#'
vis.multiple.pt <- function(x, # use pt objects directly? instead of plots? and plot internally
                            to,
                            labels,
                            legend.title) {
  
  # For cmd check
  . <- time <- state <- prob <- CI_low <- CI_upp <- NULL
  
  # Check x are ggplots
  cond <- sapply(x, function(p) "ggplot" %in% class(p))
  if (!all(cond)) stop("x should be a list of ggplots")
  
  # Check labels
  if (missing(labels)) labels <- paste0("pt_", 1:length(x))
  if (missing(legend.title)) legend.title <- "PT"
  
  # Add labels
  dfs_labelled <- lapply(1:length(x), function(i) {
    df <- x[[i]]$data
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
    ggplot2::geom_line(ggplot2::aes(col = PT), size = 2) +
    ggplot2::geom_ribbon(alpha = 0.5, fill = "lightgray") +
    ggplot2::guides(colour = ggplot2::guide_legend(title = legend.title))
  
  return(p)
}