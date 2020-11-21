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
#'         Chavent (2007). It is `1 - (sum(current inertia)/inertial[1])`.}
#'       \item{alt}{A nested tibble of alternate splits at a node. It contains
#'         `bipartsplitrow` and `bipartsplitcol` with the same meaning above.
#'         Note that this is only for information purpose. Currently `monoClust`
#'         does not support choosing an alternate splitting route. Running
#'         [MonoClust()] with `nclusters = 2` step-by-step can be run if
#'         needed.}
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
#'   \item{alt}{Indicator of having an alternate splitting route occurred when
#'     splitting.}
#'   \item{circularroot}{List of values designed for circular variable in the
#'     data set. `var` is the name of circular variable and `cut` is its first
#'     best split value. If circular variable is not available, both objects are
#'     NULL.}
#' }
#' @references
#' * Chavent, M., Lechevallier, Y., & Briant, O. (2007). DIVCLUS-T: A monothetic
#' divisive hierarchical clustering method. Computational Statistics & Data
#' Analysis, 52(2), 687-701. \doi{10.1016/j.csda.2007.03.013}.
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
                     alt = list(
                       tibble::tibble(bipartsplitrow = numeric(),
                                      bipartsplitcol = numeric()))) {

  one_row_table <- tibble::tibble(
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
