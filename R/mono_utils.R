#' Monothetic Clustering Trees Object
#'
#' The structure and objects contained in MonoClust, an object returned from
#' the [MonoClust()] function and used as the input in other functions in the
#' package.
#'
#' @name MonoClust.object
#'
#' @return
#' \describe{
#'   \item{frame}{Data frame in the form of a [tibble::tibble()] representing
#'     a tree structure with one row for each node. The columns include:
#'     \describe{
#'       \item{number}{Index of the node. Depth of a node can be derived by
#'         `number %/% 2`.}
#'       \item{var}{Name of the variable used in the split at a node or
#'         `"<leaf>"` if it is a leaf node.}
#'       \item{cut}{Splitting value, so values of `var` that are smaller than
#'         that go to left branch while values greater than that go to the right
#'         branch.}
#'       \item{n}{Cluster size, the number of observations in that cluster.}
#'       \item{inertia}{Inertia value of the cluster at that node.}
#'       \item{bipartsplitrow}{Position of the next split row in the data set
#'         (that position will belong to left node (smaller)).}
#'       \item{bipartsplitcol}{Position of the next split variable in the data
#'         set.}
#'       \item{inertiadel}{Proportion of inertia value of the cluster at that
#'         node to the inertia of the root.}
#'       \item{medoid}{Position of the data point regarded as the medoid of
#'         its cluster.}
#'       \item{loc}{y-coordinate of the splitting node to facilitate showing
#'         on the tree. See [plot.MonoClust()] for details.}
#'       \item{split.order}{Order of the splits with root is 0.}
#'       \item{inertia_explained}{Percent inertia explained as described in
#'         Chavent (2007).}
#'       \item{alt}{Indicator of an alternative cut yielding the same reduction
#'         in inertia at that split.}
#'     }}
#'   \item{membership}{Vector of the same length as the number of rows in the
#'     data, containing the value of `frame$number` corresponding to the leaf
#'     node that an observation falls into.}
#'   \item{dist}{Distance matrix calculated using the method indicated in
#'     `distmethod` argument of [MonoClust()].}
#'   \item{terms}{Vector of variable names in the data that were used to split.}
#'   \item{centroids}{Data frame with one row for centroid value of each
#'     cluster.}
#'   \item{medoids}{Named vector of positions of the data points regarded as
#'     medoids of clusters.}
#'   \item{circularroot}{List of values designed for circular variable in the
#'     data set. `var` is the name of circular variable and `cut` is its first
#'     best split value. If circular variable is not available, both objects are
#'     NULL.}
#' }
#' @references
#' * Chavent, M., Lechevallier, Y., & Briant, O. (2007). DIVCLUS-T: A monothetic
#' divisive hierarchical clustering method. Computational Statistics & Data
#' Analysis, 52(2), 687–701. <doi:10.1016/j.csda.2007.03.013>.
#' @seealso [MonoClust()].
NULL

#' Find the Closest Cut
#'
#' Find the cuts for a quantitative variable. These cuts are what we are
#' going to consider when thinking about bi-partitioning the data. For a
#' quantitative column, find the next larger value of each value, if it is the
#' largest, that value + 1
#'
#' @param col a quantitative vector.
#'
#' @return a quantitative vector which contains the closest higher cut.
#' @keywords internal
find_closest <- function(col) {
  purrr::map_dbl(col, ~ ifelse(.x == max(col),
                               .x + 1,
                               min(col[which(col - .x > 0)])))
}

#' Add/Subtract Circular Values in Degrees/Radian
#'
#' Add/subtract two circular variables in degrees (`%cd+%` and `%cd-%`) and
#' radian (`%cr+%` and `%cr-%`).
#'
#' @param x,y Circular values in degrees/radians.
#'
#' @return A value between [0, 360) in degrees or [0, 2*pi) in radian.
#' @name circ_arith
#' @examples
#' 90 %cd+% 90
#'
#' 250 %cd+% 200
#'
#' 25 %cd-% 80
#'
#' pi %cr+% (pi/2)
#'
NULL

#' @export
#' @rdname circ_arith
`%cd+%` <- function(x, y) {
  return((x + y) %% 360)
}

#' @export
#' @rdname circ_arith
`%cd-%` <- function(x, y) {
  return((x - y) %% 360)
}

#' @export
#' @rdname circ_arith
`%cr+%` <- function(x, y) {
  return((x + y) %% (2 * pi))
}

#' @export
#' @rdname circ_arith
`%cr-%` <- function(x, y) {
  return((x - y) %% (2 * pi))
}

#' Cluster Inertia Calculation
#'
#' Calculate inertia for a given subset of the distance matrix from the original
#' data set provided to X. Assumes that distance matrices are stored as matrices
#' and not distance objects
#'
#' @param X distance matrix, not an object of some distance measure
#'
#' @return inertia value of the matrix, formula in Chavent (1998). If X is a
#'   single number, return 0.
#' @keywords internal
inertia_calc <- function(X) {

  if (!is.numeric(X) && !is.matrix(X))
    stop("X has to be a numerical value or matrix.")

  # If singleton cluster, inertia is 0
  inertia_value <- ifelse(length(X) == 1 && is.numeric(X),
                          0,
                          sum(X^2) / (dim(X)[1] * 2))
  return(inertia_value)
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
#' @keywords internal
medoid <- function(members, dist_mat) {
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

#' Create A New Node for Split Data Frame
#'
#' This function is just a helper to make sure that the default values of the
#' split data frame is correct when unspecified. It helps reduce type error,
#' especially when moving to use dplyr which is stricter in data types.
#'
#' @param number Row index of the data frame.
#' @param var Whether it is a leaf, or the name of the next split variable.
#' @param cut The splitting value, so values (of `var`) smaller than that
#'   go to left branch while values greater than that go to right branch.
#' @param n Cluster size. Number of observations in that cluster.
#' @param inertia Inertia value of the cluster at that node.
#' @param bipartsplitrow Position of the next split row in the data set (that
#'   position will belong to left node (smaller)).
#' @param bipartsplitcol Position of the next split variable in the data set.
#' @param inertiadel The proportion of inertia value of the cluster at that node
#'   to the inertia of the root.
#' @param medoid Position of the data point regarded as the medoid of its
#'   cluster.
#' @param loc y-coordinate of the splitting node to facilitate showing on the
#'   tree. See [plot.MonoClust()] for details.
#' @param split.order Order of the splits. Root is 0, and increasing.
#' @param inertia_explained Percent inertia explained as described in Chavent
#'   (2007)
#' @param alt Indicator of an alternative cut yielding the same reduction in
#'   inertia at that split.
#'
#' @references
#' * Chavent, M., Lechevallier, Y., & Briant, O. (2007). DIVCLUS-T: A monothetic
#' divisive hierarchical clustering method. Computational Statistics & Data
#' Analysis, 52(2), 687–701. https://doi.org/10.1016/j.csda.2007.03.013
#'
#' @return A tibble with only one row and correct default data type for even an
#'   unspecified variables.
#' @keywords internal
new_node <- function(number,
                     var,
                     cut = -99L,
                     n,
                     inertia,
                     bipartsplitrow = -99L,
                     bipartsplitcol = -99L,
                     inertiadel = 0,
                     inertia_explained = -99,
                     medoid,
                     loc,
                     split.order = -99L,
                     alt = FALSE) {

  one_row_table <- dplyr::tibble(
    number, var, cut, n,
    inertia, bipartsplitrow,
    bipartsplitcol, inertiadel,
    inertia_explained, medoid,
    loc,
    split.order,
    alt)

  return(one_row_table)
}

#' Find Centroid of the Cluster
#'
#' Centroid is point whose coordinates are the means of their cluster values.
#'
#' @inheritParams checkem
#'
#' @return A data frame with coordinates of centroids
#' @keywords internal
centroid <- function(data, frame, cloc) {

  leaves <- frame$number[frame$var == "<leaf>"]
  names(leaves) <- leaves
  centroid_list <- vector("list", length(leaves))

  centroid_list <-
    purrr::map_dfr(leaves, function(x) {
      cluster <- data[cloc == x, ]
      centroid <- purrr::map_dbl(cluster, mean)
      return(centroid)
    },
    .id = "cname")

  return(centroid_list)
}

#' Find Tree Depth Based on Node Indices
#'
#' @param nodes Vector of node indices in the tree.
#'
#' @details
#' When building MonoClust tree, the node index was created with the rule that
#'   new node indices are the split node times 2 plus 0 (left) and 1 (right).
#'   Therefore, this function is just a back-transform, taking a log base 2.
#'
#' @return Depth of the node, with 0 is the root relative to the input.
#'
#' @keywords internal
tree_depth <- function(nodes) {
  depth <- floor(log2(nodes) + 1e-07)
  return(depth - min(depth))
}
