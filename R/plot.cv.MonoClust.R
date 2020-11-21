#' Plot the Mean Square Error with Error Bar for +/- 1 Standard Error
#'
#' @param x A `cv.MonoClust` object (output of [cv.test()]).
#' @param main Overall title for the plot.
#' @param xlab Title for x axis.
#' @param ylab Title for y axis.
#' @param type What type of plot should be drawn. See [graphics::par()].
#' @param lty The line type.
#' @param err.col Color of the error bars.
#' @param ... Arguments to be passed to [graphics::plot.default()].
#' @param err.width Width of the bars.
#'
#' @return A line plot with error bars.
#' @seealso Plot using ggplot2 [ggcv()]
#' @export
#'
#' @examples
#' \donttest{
#' library(cluster)
#' data(ruspini)
#'
#' # 10-fold cross-validation
#' cptable <- cv.test(ruspini, minnodes = 2, maxnodes = 4)
#' plot(cptable)
#' }
plot.cv.MonoClust <- function(x,
                              main = "MSE for CV of monothetic clustering",
                              xlab = "Number of clusters",
                              ylab = "MSE +/- 1 SE",
                              type = "b",
                              lty = 2,
                              err.col = "red",
                              err.width = 0.1,
                              ...) {

  if (missing(x))
    stop("\"x\" is required.")
  if (!inherits(x, "cv.MonoClust"))
    stop("Not a legitimate \"cv.MonoClust\" object")

  graphics::plot(x$cv$ncluster, x$cv$MSE, type = type, lty = lty,
                 ylim = c(0, max(x$cv$MSE + x$cv$`Std. Dev.`)),
                 main = main, ylab = ylab, xlab = xlab, ...)

  error_bar(x$cv$ncluster, x$cv$MSE, x$cv$`Std. Dev.`, col = err.col)

  return(invisible(x))
}

#' Make Error Bars
#'
#' @param x x coordinates.
#' @param y y coordinates.
#' @param upper Distance from y to the upper bar.
#' @param lower Distance from y to the lower bar.
#' @param length Length of the horizontal bar.
#' @param ... Other arguments to [graphics::arrows()]
#'
#' @return Plot
#'
#' @keywords internal
error_bar <- function(x, y, upper, lower = upper, length = 0.1, ...) {
  if (length(x) != length(y) |
      length(y) != length(lower) |
      length(lower) != length(upper)) stop("vectors must be same length")
  graphics::arrows(x, y + upper,
                   x, y - lower,
                   angle = 90, code = 3, length = length, ...)
  return(invisible(x))
}
