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
predict.MonoClust <- function(object, newdata, type = c("centroid", "medoid")) {

  ####### DELETE: For debugging purpose ########### require(cluster) data(ruspini) source('MonoClust.R')
  ####### object <- MonoClust(ruspini, nclusters = 4) newdata <- list(x = c(48, 90), y = c(93, 70))
  ####### na.action <- na.pass

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

    # node <- frame$number
    # depth <- tree_depth(node)

    # It would be better to create a jump table for reference of tree walking
    jump_table <- dplyr::select(frame, number, bipartvar, cut)
    jump_table <- tibble::add_column(jump_table, left = NA, right = NA)

    if (nrow(jump_table) >= 2) {
      jump_table$left <- match(jump_table$number * 2, jump_table$number)
      jump_table$right <- match(jump_table$number * 2 + 1, jump_table$number)
    } else jump_table[, c("bipartvar", "number")] <- NA

    # Now tracing the tree to find the cluster
    terminal_node <- apply(new_data, 1, tree_walk, jump_table = jump_table)
  }
  # ylevels <- attr(object, "ylevels")
  # nclass <- length(ylevels)
  # if (missing(type) && nclass > 0L)
  #   type <- "prob"
  if (type == "centroid") {
    centroids <- object$centroids
    ret <- centroids[match(terminal_node, centroids$cname), ]
  } else {
    ret <- tibble::tibble(cname = terminal_node,
                  medoid = object$frame$medoid[terminal_node])
  }
  # } else if (type == 'matrix') { pred <- frame$yval2[terminal_node, ] dimnames(pred) <-
  # list(names(terminal_node), NULL) } else if (type == 'class' && nclass > 0L) { if (length(ylevels) ==
  # 0L) stop('type 'class' is only appropriate for classification') pred <-
  # factor(ylevels[frame$yval[terminal_node]], levels = ylevels) names(pred) <- names(terminal_node) } else if
  # (type == 'prob' && nclass > 0L) { pred <- frame$yval2[terminal_node, 1L + nclass + 1L:nclass, drop =
  # FALSE] dimnames(pred) <- list(names(terminal_node), ylevels) } else stop('Invalid prediction for
  # \'rpart\' object')

  # Expand out the missing values in the result But only if operating on the original dataset if
  # (missing(newdata) && !is.null(object$na.action)) pred <- naresid(object$na.action, pred)
  return(ret)
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
    var <- jump_table$bipartvar[current_node]
    value <- new_point[var]

    current_node <- dplyr::if_else(value < jump_table$cut[current_node],
                                   jump_table$left[current_node],
                                   jump_table$right[current_node])
  }

  # Extract the leaf info out. TODO: Output the necessary info.
  # QUESTION: Do we need to re-calculate any properties? Like medoid?
  node <- jump_table$number[current_node]
  names(node) <- current_node
  return(node)
}
