#' Cluster Inertia Calculation
#'
#' Calculate inertia for a given subset of the distance matrix from the original
#' data set provided to `x`. Assumes that distance matrices are stored as
#' matrices and not distance objects.
#'
#' @param x Distance matrix, not an object of some distance measure.
#'
#' @return Inertia value of the matrix, formula in Chavent (1998). If `x` is a
#'   single number, return 0.
#' @export
inertia_calc <- function(x) {

  if (!is.numeric(x) && !is.matrix(x))
    stop("x has to be a numerical value or matrix.")

  # If singleton cluster, inertia is 0
  inertia_value <- ifelse(length(x) == 1 && is.numeric(x),
                          0,
                          sum(x^2) / (dim(x)[1] * 2))
  return(inertia_value)
}
