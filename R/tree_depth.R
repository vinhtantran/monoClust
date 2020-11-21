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
