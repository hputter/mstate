# Helper function to set default plot colours
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
    
  } else if(type == "areas") { # areas
    # Viridis colours, no limit
    cols <- viridis::viridis_pal()(n)
  } else stop("Type should be either 'lines' or 'areas'!")
  
  return(cols)
}