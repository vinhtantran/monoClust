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
