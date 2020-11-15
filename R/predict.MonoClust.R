#' Predictions from a MonoClust Object
#'
#' Predict the cluster memberships of a new data set from a `MonoClust` object.
#'
#' @param object MonoClust result object.
#' @param newdata Data frame containing the values to be predicted. If missing,
#'   the memberships of the MonoClust object are returned.
#' @param type Type of returned cluster representatives. Either `"centroid"` to
#'   return the centroid values of the terminal clusters, or `"medoid"` to
#'   return the index of the medoid observations in the clustered data set.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A tibble of cluster index in `cname` and either centroid values or
#'   medoid observations index based on the value of `type` argument.
#'
#' @seealso [predict]
#' @export
#'
#' @importFrom stats na.pass
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#'
#' set.seed(1234)
#' test_index <- sample(1:nrow(ruspini), nrow(ruspini)/5)
#' train_index <- setdiff(1:nrow(ruspini), test_index)
#' ruspini_train <- ruspini[train_index, ]
#' ruspini_test <- ruspini[test_index, ]
#'
#' ruspini_train_4sol <- MonoClust(ruspini_train, nclusters = 4)
#' predict(ruspini_train_4sol, newdata = ruspini_test)
predict.MonoClust <- function(object, newdata, type = c("centroid", "medoid"),
                              ...) {

  if (!inherits(object, "MonoClust"))
    stop("Not a legitimate \"MonoClust\" object")

  type <- match.arg(type)

  terminal_node <- NA
  if (missing(newdata)) {
    terminal_node <- object$Membership
  } else {
    frame <- object$frame
    new_data <- tibble::tibble(newdata)
    # Check to see if the newdata has the same structure as the object
    terms <- object$terms
    new_terms <- colnames(new_data)
    if (!length(intersect(terms, new_terms)) == length(terms))
      stop("The new data set does not have the same variables as the original
           data set")

    # It would be better to create a jump table for reference of tree walking
    jump_table <- make_jump_table(frame)

    # Now tracing the tree to find the cluster
    terminal_node <- apply(new_data, 1, tree_walk, jump_table = jump_table)
  }

  if (type == "centroid") {
    centroids <- object$centroids
    ret <- centroids[match(terminal_node, centroids$cname), ]
  } else {
    ret <- tibble::tibble(cname = terminal_node,
                          medoid = object$frame$medoid[terminal_node])
  }

  return(ret)
}

#' Create Jump Table
#'
#' Create jump table from the MonoClust's frame object. `number` and `var` will
#' be used to create the table.
#'
#' @param frame MonoClust's frame object
#'
#' @return Jump table with `number`, `var`, and two new columns `left` and
#'   `right` indicate the left and right number at split.
#' @keywords internal
make_jump_table <- function(frame) {
  jump_table <- frame[, c("number", "var", "cut")]
  jump_table <- tibble::add_column(jump_table, left = NA, right = NA)

  if (nrow(jump_table) >= 2) {
    jump_table$left <- match(jump_table$number * 2, jump_table$number)
    jump_table$right <- match(jump_table$number * 2 + 1, jump_table$number)
  } else
    # Special case when tree didn't split, ncluster = 1
    jump_table[, "var"] <- NA

  return(jump_table)
}

#' Traverse a Tree to Find the Leaves (Terminal Nodes)
#'
#' @param new_point New data point
#' @param jump_table Jump table
#'
#' @return The index of the terminal node after traversing the new data point on
#'   the tree.
#'
#' @keywords internal
tree_walk <- function(new_point, jump_table) {
  current_node <- 1
  while (!is.na(jump_table$cut[current_node])) {
    # If it's not a leaf node, trace
    var <- jump_table$var[current_node]
    value <- new_point[var]

    current_node <- dplyr::if_else(value < jump_table$cut[current_node],
                                   jump_table$left[current_node],
                                   jump_table$right[current_node])
  }

  # Extract the leaf info out.
  node <- jump_table$number[current_node]
  names(node) <- current_node
  return(node)
}
