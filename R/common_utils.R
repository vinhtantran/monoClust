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

#' Find Tree Depth Based on Node Indexes
#'
#' @param nodes Vector of node indexes in the tree.
#'
#' @details
#' When building MonoClust tree, the node index was created with the rule that
#'   new node indexes are the split node times 2 plus 0 (left) and 1 (right).
#'   Therefore, this function is just a back-transform, taking a log base 2.
#'
#' @return Depth of the node, with 0 is the root relative to the input.
#'
#' @export
tree_depth <- function(nodes) {
  if (!is.numeric(nodes))
    stop("\"node\" has to be a numerical value.")

  depth <- floor(log2(nodes) + 1e-07)
  return(depth - min(depth))
}

#' Find Medoid of the Cluster
#'
#' Medoid is the point that has minimum distance to all other points in the
#' cluster.
#'
#' @param members index vector indicating which observation belongs to the
#'   cluster.
#' @param dist_mat distance matrix of the whole data set.
#'
#' @return index of the medoid point in the members vector.
#' @export
medoid <- function(members, dist_mat) {
  if (max(members) > dim(dist_mat)[1])
    stop("Value(s) of \"members\" is out of bound of \"dist_mat\".")

  index <- NULL

  if (length(members) == 0) {
    index <- 0
  } else if (length(members) == 1) {
    index <- members
  } else {
    dists <- purrr::map_dbl(purrr::array_branch(dist_mat[members, members], 1),
                            sum)
    medoid <- members[which(dists == min(dists))]
    index <- medoid[1]
  }

  return(index)
}

#' Coerce Similar Object to MonoClust
#'
#' The function turns a MonoClust-similar object into MonoClust object so it
#' can use supported functions for MonoClust such as [print.MonoClust()] and
#' [plot.MonoClust()].
#'
#' `as_MonoClust()` is an S3 generic. The function itself doesn't run unless
#' it is implemented for another similar object. Currently, this function is not
#' implemented within `monoClust` package.
#'
#' @param x An object that can be coerced to MonoClust object.
#' @param ... For extensibility.
#'
#' @export
as_MonoClust <- function(x, ...) {
  UseMethod("as_MonoClust")
}

#' @export
#' @rdname as_MonoClust
as_MonoClust.default <- function(x, ...) {
  warning(paste("as_MonoClust does not know how to handle object of class ",
                class(x), "."))
}
