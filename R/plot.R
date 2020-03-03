#' Plots the summary matrix
#' 
#' Plots the summary matrix with counts on y-axis and classification on x-axis.
#' 
#' @param x Summary matrix returned from dbcompare
#' @param log Specifies whether log(Counts) should be plotted (default)
#' @param las Direction of the labels on x-axis. Default is 3 which gives
#' perpendicular labels
#' @param xlab Axis label
#' @param ylab Axis label
#' @param ... Other plot options
#' @return A plot of the summary matrix. The counts are on log10 scale and the
#' x-axis is labeled by appropriate matching/partially-matching levels.
#' @author James Curran and Torben Tvedebrink
#' @seealso dbCompare,print.dbcompare
#' @examples
#' 
#'   \dontrun{
#'   data(dbExample)
#'   M = dbCompare(dbExample,hit=5)
#'   plot(M)
#'   }
#' 
#' @rawNamespace S3method(plot, dbcompare)
#' @importFrom graphics plot axis box
#' @export plot.dbcompare
plot.dbcompare <- function(x, log = "y", las = 3, xlab = "Match/Partial", ylab = "Counts", ...) {
  nl <- attributes(x)$call$loci
  levs <- dbCats(nl, vector = TRUE)
  if (is.matrix(x$m)) 
    mvec <- t(x$m)[up.tri(x$m)] else mvec <- x$m
  mvec[mvec == 0] <- NA
  if (attributes(x)$call$collapse) {
    mcol <- ifelse(x$m == 0, NA, x$m)
    if (xlab == "Match/Partial") 
      xlab <- "Total number of matching alleles"
    graphics::plot(0:(length(mcol) - 1), mcol, log = log, xlab = xlab, ylab = ylab, ...)
  } else {
    graphics::plot(1:length(levs), mvec, axes = FALSE, xlab = xlab, ylab = ylab, log = log, ...)
    graphics::axis(1, at = 1:length(levs), labels = levs, las = las)
    graphics::axis(2)
    graphics::box()
  }
}
