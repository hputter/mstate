# Helper function to set default plot colours
set_colours <- function(n, type) {
  
  if (type == "lines") {
    
    # Rcolorbrewer Dark2 - max 8 colours
    full_pal <- RColorBrewer::brewer.pal(8, "Dark2")
    
    if (n <= 8) {
      cols <- full_pal[1:n]
    } else {
      dark2_extended <- grDevices::colorRampPalette(full_pal)
      cols <- dark2_extended(n)
    }
    
  } else if(type == "areas") { 
    cols <- viridisLite::viridis(n = n)
  } else stop("Type should be either 'lines' or 'areas'!")
  
  return(cols)
}