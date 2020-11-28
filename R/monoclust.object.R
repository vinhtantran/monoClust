#' Monothetic Clustering Tree Object
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
