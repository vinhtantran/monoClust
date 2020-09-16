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

#' Add/Subtract From A Circular Value
#'
#' Add a value to a circular value or variable. When the value reaches 360 or 0,
#' it will become 0.
#'
#' @param variable circular variable (in degree 0-360).
#' @param shift the added value to the circular variable, can be positive or
#'   negative if want to subtract.
#'
#' @return shifted circular value or variable.
#' @keywords internal
cshift <- function(variable, shift) {
  # In case shift is larger a full circle
  shift_less_360 <- shift %/% 360
  # Add value
  variable_shifted <- variable + shift_less_360
  variable_shifted <- ifelse(variable_shifted < 0,
                             variable_shifted + 360,
                             variable_shifted)
  variable_shifted <- ifelse(variable_shifted >= 360,
                             variable_shifted - 360,
                             variable_shifted)
  return(variable_shifted)
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

#' Circular Distance using Gower's
#'
#' calculates the distance matrix within a circular variable using Gower's
#' distance. Written by Garland Will.
#'
#' @param x a numeric vector of circular values
#'
#' @return object of class "dist"
#' @keywords internal
#' @importFrom stats as.dist
circ_dist <- function(x) {
  # Assumes x is just a single variable
  dist1 <- matrix(0, nrow = length(x), ncol = length(x))
  for (i in seq_len((length(x) - 1))) {
    for (j in (i+1):length(x)) {
      dist1[j, i] = min(abs(x[i] - x[j]), (360 - abs(x[i] - x[j])))/180
    }
  }
  return(stats::as.dist(dist1))
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
                     split.order = -99L) {

  one_row_table <- dplyr::tibble(
    number, var, cut, n,
    # wt,
    inertia, bipartvar, bipartsplitrow,
    bipartsplitcol, inertiadel, yval, inertia_explained, medoid,
    # category,
    loc,
    split.order)

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
