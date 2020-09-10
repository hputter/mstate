#' @method summary Cuminc
#' @export
summary.Cuminc <- function(object, ...)
{
  if (!inherits(object, "Cuminc"))
    stop("'object' must be a 'Cuminc' object")
  
  summary(attr(object, "survfit"), ...)
}
