#' Summary method for a summary.Cuminc object
#' 
#' @param object Object of class 'Cuminc', to be summarised
#' @param \dots Further arguments to summarise
#'
#' @export 
summary.Cuminc <- function(object, ...)
{
  if (!inherits(object, "Cuminc"))
    stop("'object' must be a 'Cuminc' object")
  
  summary(attr(object, "survfit"), ...)
}

