#' @import survival
#' @importFrom data.table .N `:=` .SD melt.data.table data.table setDT
#' @importFrom magrittr `%>%`
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
#' @importFrom graphics box lines par plot polygon text title
#' @importFrom stats as.formula delete.response model.frame model.matrix 
#' model.offset model.response rnorm terms time aggregate approxfun 
#' quantile rpois weighted.mean approx qnorm
#' @importFrom utils flush.console head tail
#' @importFrom lattice xyplot
#' 
#' @useDynLib mstate, .registration=TRUE
NULL