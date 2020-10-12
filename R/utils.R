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

#' Add/Subtract Circular Values in Degrees
#'
#' Add/subtract two circular variables.
#'
#' @param x,y Circular values in degrees.
#'
#' @return A value between [0, 360).
#' @name circ_arith
#' @examples
#' 90 %circ+% 90
#'
#' 250 %circ+% 200
#'
#' 25 %circ-% 80
#'
NULL

#' @export
#' @rdname circ_arith
`%circ+%` <- function(x, y) {
  return((x + y) %% 360)
}

#' @export
#' @rdname circ_arith
`%circ-%` <- function(x, y) {
  return((x - y) %% 360)
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
  # there are cases when a cluster has only 1 point, say, 1st point, then
  # dist[1,1] is a numeric value, not matrix.
  #MG, 9/25: Should this then return a value of 0 for inertia? If you go back to
  # (y-mean(y))^2, then maybe set the return to 0?
  if (!is.numeric(X) && !is.matrix(X))
    stop("X has to be a numerical value or matrix.")

  # If singleton cluster, inertia is 0
  inertia_value <- ifelse(length(X) == 1 && is.numeric(X),
                          0,
                          sum(X^2) / (dim(X)[1] * 2))
  return(inertia_value)
}

#' Distance Matrix of Circular Variables
#'
#' Calculates the distance matrix of observations with circular variables using
#' an adapted version of Gower's distance. This distance should be compatible
#' with the Gower's distance for other variable types.
#'
#' @param frame A data frame with all columns are circular.
#'
#' @details
#' The distance between two observations {i} and {j} of a circular variable {q}
#' is suggested to be
#'
#' \deqn{(y_{iq}, y_{jq}) = \frac{180 - |180 - |y_{iq} - y_{jq}||}{180}.}
#'
#' @return Object of class "dist".
#'
#' @importFrom stats as.dist
#'
#' @references
#' * Tran, T. V. (2019). Chapter 3. Monothetic Cluster Analysis with Extensions
#' to Circular and Functional Data. Montana State University - Bozeman.
#' @export
circ_dist <- function(frame) {
  # Assumes x is a data frame with columns are all circular variables
  # TODO Extend it to more than one circular variable

  if (is.null(frame))
    stop("frame has to be a data set with all columns are circular.")

  frame <- frame %circ+% 0

  gower_circ <- function(x, y) abs(180 - abs(180 - abs(x - y))) / 180

  list_dist <-
    purrr::map(frame, function(x) {

      dist <- matrix(0, ncol = length(x), nrow = length(x))

      for (i in seq_len(length(x) - 1))
        for (j in (i + 1):length(x))
          dist[j, i] <- gower_circ(x[i], x[j])

      return(dist)


      # pairs_of_obs <- purrr::cross2(x, x)
      #
      # dist_flat <- purrr::flatten_dbl(
      #   purrr::map(pairs_of_obs,
      #              ~ dplyr::if_else(.x[[1]] > .x[[2]],
      #                               gower_circ(.x[[1]], .x[[2]]),
      #                               0)))
      # return(dist_flat)
    })

  ret <- matrix(purrr::pmap_dbl(list_dist, sum)/length(list_dist), ncol = nrow(frame))
  return(as.dist(ret))
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
#' @param cut The splitting value, so values (of `bipartvar`) smaller than that
#'   go to left branch while values greater than that go to right branch.
#' @param n Cluster size. Number of observations in that cluster.
#' @param inertia Inertia value of the cluster at that node.
#' @param bipartvar Name of the next split variable, match `var` if `var` is not
#'   a leaf.
#' @param bipartsplitrow Position of the next split row in the data set (that
#'   position will belong to left node (smaller)).
#' @param bipartsplitcol Position of the next split variable in the data set.
#' @param inertiadel The proportion of inertia value of the cluster at that node
#'   to the inertia of the root.
#' @param yval In rpart, it's the estimated response value. Here in MonoClust,
#'   it has no meaning. Will be removed in the future.
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
                     # Remove it is not implemented
                     # wt,
                     inertia,
                     bipartvar = "",
                     bipartsplitrow = -99L,
                     bipartsplitcol = -99L,
                     inertiadel = 0,
                     # TODO: Replace yval by inertia_explained
                     yval,
                     inertia_explained = -99,
                     medoid,
                     # Remove later because there's no categorical
                     # category = 0,
                     loc,
                     split.order = -99L,
                     alt = FALSE) {

  one_row_table <- dplyr::tibble(
    number, var, cut, n,
    # wt,
    inertia, bipartvar, bipartsplitrow,
    bipartsplitcol, inertiadel, yval, inertia_explained, medoid,
    # category,
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

  # MODIFY: Tan, 9/9/20. Remove categorical variable for now.
  ## ADD, Tan, 12/15, function to calculate the mean of each cluster.
  ## Currently do not work for categorical variables
  # Don't calculate if there is qualitative variable
  # if (qualtog) NA

  leaves <- frame$number[frame$var == "<leaf>"]
  names(leaves) <- leaves
  centroid.list <- vector("list", length(leaves))

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

#' What to Use with ForEach
#'
#' @param x A binary output of [getDoParWorkers()].
#'
#' @return Appropriate operator depending on whether parallel processing is
#'   activated or not.
#' @importFrom foreach `%dopar%`
#' @importFrom foreach `%do%`
#' @keywords internal
getOper <- function(x) {
  if(x) `%dopar%` else `%do%`
}
