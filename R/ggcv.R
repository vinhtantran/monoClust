#' GGPlot the Mean Square Error with Error Bar for +/- 1 Standard Error
#'
#' @param cv.obj A `cv.MonoClust` object (output of [cv.test()]).
#' @param title Overall title for the plot.
#' @param xlab Title for x axis.
#' @param ylab Title for y axis.
#' @param type What type of plot should be drawn. Choosing between `"l"` (line
#'   only), `"p"` (point only), and `"b"` (both line and point).
#' @param linetype The line type. See `vignette("ggplot2-specs")`.
#' @param err.col Color of the error bars.
#' @param err.width Width of the bars.
#'
#' @import dplyr
#' @importFrom rlang .data
#' @seealso Plot using base R [plot.cv.MonoClust()]
#' @export
#'
#' @examples
#' \donttest{
#' library(cluster)
#' data(ruspini)
#'
#' # 10-fold cross-validation
#' cptable <- cv.test(ruspini, minnodes = 2, maxnodes = 4)
#' ggcv(cptable)
#' }
ggcv <- function(cv.obj,
                 title = "MSE for CV of monothetic clustering",
                 xlab = "Number of clusters",
                 ylab = "MSE +/- 1 SE",
                 type = c("b", "p", "l"),
                 linetype = 2,
                 err.col = "red",
                 err.width = 0.2) {

  if (missing(cv.obj))
    stop("\"cv.obj\" is required.")
  if (!inherits(cv.obj, "cv.MonoClust"))
    stop("Not a legitimate \"cv.MonoClust\" object")
  type <- match.arg(type)

  cv_table <- cv.obj$cv

  p <-
    cv_table %>%
    mutate(upper1SD = .data$MSE + .data$`Std. Dev.`,
           lower1SD = .data$MSE - .data$`Std. Dev.`) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$ncluster, y = .data$MSE)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower1SD,
                                        ymax = .data$upper1SD),
                           color = err.col, width = err.width) +
    ggplot2::scale_x_continuous(breaks = seq_len(nrow(cv_table))) +
    ggplot2::labs(x = xlab,
                  y = ylab,
                  title = title)
  if (type %in% c("p", "b"))
    p <- p + ggplot2::geom_point()
  if (type %in% c("l", "b"))
    p <- p + ggplot2::geom_line(linetype = linetype)

  return(p)
}
