#' Distance Matrix of Circular Variables
#'
#' Calculates the distance matrix of observations with circular variables using
#' an adapted version of Gower's distance. This distance should be compatible
#' with the Gower's distance for other variable types.
#'
#' @param frame A data frame with all columns are circular measured in degrees.
#'
#' @details
#' The distance between two observations {i} and {j} of a circular variable {q}
#' is suggested to be
#'
#' \deqn{(y_{iq}, y_{jq}) = \frac{180 - |180 - |y_{iq} - y_{jq}||}{180}.}
#'
#' @return Object of class "dist".
#'
#' @seealso [stats::dist()]
#'
#' @references
#' * Tran, T. V. (2019). Chapter 3. Monothetic Cluster Analysis with Extensions
#' to Circular and Functional Data. Montana State University - Bozeman.
#' @export
#' @examples
#' # Make a sample data set of 20 observations with 2 circular variables
#' data <- data.frame(var1 = sample.int(359, 20),
#'                    var2 = sample.int(359, 20))
#' circ_dist(data)
circ_dist <- function(frame) {

  if (missing(frame))
    stop("frame has to be a data set with all columns are circular.")

  frame <- frame %cd+% 0

  gower_circ <- function(x, y) abs(180 - abs(180 - abs(x - y))) / 180

  list_dist <-
    purrr::map(frame, function(x) {

      dist <- matrix(0, ncol = length(x), nrow = length(x))

      for (i in seq_len(length(x) - 1))
        for (j in (i + 1):length(x))
          dist[j, i] <- gower_circ(x[i], x[j])

      return(dist)
    })

  ret <- matrix(purrr::pmap_dbl(list_dist, sum) / length(list_dist),
                ncol = nrow(frame))
  return(stats::as.dist(ret))
}
