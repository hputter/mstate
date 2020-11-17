#' Automatic choice of default colours for plot
#'
#' @param n Number of groups/lines/areas to be plotted. Generally number
#' of states or transitions 
#' @param type Either "line" or "areas"
#'
#' @return Vector of colours
set_colours <- function(n, type) {
  
  if (type == "lines") {
    
    # Rcolorbrewer Dark2 - max 8 colours
    if (n <= 8) {
      full_pal <- RColorBrewer::brewer.pal(8, "Dark2")
      cols <- full_pal[1:n]
    } else {
      dark2_extended <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Dark2")
      )
      cols <- dark2_extended(n)
    }
    
  } else {
    # Viridis colours, no limit
    cols <- viridis::viridis_pal()(n)
  }
  
  return(cols)
}