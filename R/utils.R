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
  if (!is.numeric(X) && !is.matrix(X)) stop("X has to be a numerical value or matrix.")

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
circ_dist <- function(x) {
  # Assumes x is just a single variable
  dist1 <- matrix(0, nrow = length(x), ncol = length(x))
  for (i in seq_len((length(x) - 1))) {
    for (j in (i+1):length(x)) {
      dist1[j, i] = min(abs(x[i] - x[j]), (360 - abs(x[i] - x[j])))/180
    }
  }
  return(as.dist(dist1))
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
