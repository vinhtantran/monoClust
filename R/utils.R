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

#' Plot the monoClust Tree.
#'
#' This function plots the MonoClust tree. It is partially inspired by rpart
#' package.
#'
#' @inheritParams plot.MonoClust
#'
#' @return Plot of tree
#'
#' @keywords internal
plot_tree <- function(x, uniform = FALSE, branch = 1, margin = 0,
                      minbranch = 0.3, rel.loc.x = TRUE, ...) {

  # Length of the root node
  bar <- 0.03

  dev <- grDevices::dev.cur()
  if (dev == 1)
    dev <- 2

  temp <- plot_prep_node(x, uniform = uniform,
                         minbranch = minbranch)
  xx <- temp$x
  yy <- temp$y

  if (rel.loc.x)
    xx <- x$frame$loc + (abs(min(x$frame$loc))) + 1

  temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
  temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
  graphics::plot.default(temp1, temp2, type = "n", axes = FALSE, xlab = "",
                         ylab = "", ...)
  node <- x$frame$number
  temp <- plot_prep_branch(xx, yy, node, branch)
  if (branch > 0)
    graphics::lines(c(xx[1], xx[1]), c(yy[1], yy[1] + bar * diff(range(yy))),
                    ...)
  graphics::lines(c(temp$x), c(temp$y))
  invisible(list(x = xx, y = yy))
}

#' Calculate Nodes Coordinates
#'
#' @param tree MonoClust result object.
#' @inheritParams plot_tree
#'
#' @return Nodes coordinates in a list of x and y axis.
#' @keywords internal
plot_prep_node <- function(tree, uniform = FALSE, minbranch = 0.3) {
  frame <- tree$frame
  node <- as.numeric(frame$number)
  depth <- tree_depth(node)
  is_leaf <- frame$var == "<leaf>"
  if (uniform)
    y <- (1 + max(depth) - depth) / max(depth, 4)
  else {
    y <- dev <- frame$inertia
    temp <- split(seq(node), depth)

    # Index of parent of nodes in node
    parent <- match(floor(node / 2), node)

    sibling <- match(ifelse(node %% 2, node - 1, node + 1), node)

    for (i in temp[-1]) {
      temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
      y[i] <- y[parent[i]] - temp2
    }
    fudge <- minbranch * diff(range(y)) / max(depth)
    for (i in temp[-1]) {
      temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
      haskids <- !(is_leaf[i] & is_leaf[sibling[i]])
      y[i] <- y[parent[i]] - ifelse(temp2 <= fudge & haskids,
                                    fudge, temp2)
    }
    y <- y / (max(y))
  }
  x <- double(length(node))
  x[is_leaf] <- seq(sum(is_leaf))
  left_child <- match(node * 2, node)
  right_child <- match(node * 2 + 1, node)

  temp <- split(seq(node)[!is_leaf], depth[!is_leaf])
  for (i in rev(temp)) x[i] <- 0.5 * (x[left_child[i]] + x[right_child[i]])

  return(list(x = x, y = y))
}

#' Calculate Branch Coordinates
#'
#' @param x Nodes x-coordinates.
#' @param y Nodes y-coordinates.
#' @param node Nodes row number.
#' @inheritParams plot_tree
#'
#' @return Branch coordinates in a list of x and y axis.
#' @keywords internal
plot_prep_branch <- function(x, y, node, branch = 0) {
  is_left <- (node %% 2) == 0
  node_left <- node[is_left]
  parent <- match(node_left / 2, node)
  sibling <- match(node_left + 1, node)
  temp <- (x[sibling] - x[is_left]) * (1 - branch) / 2
  xx <- rbind(x[is_left], x[is_left] + temp, x[sibling] - temp,
              x[sibling], NA)
  yy <- rbind(y[is_left], y[parent], y[parent], y[sibling],
              NA)
  list(x = xx, y = yy)
}

#' Implementation of Print Labels on MonoClust Tree
#'
#' This function plots the labels onto the MonoClust tree. It is partially
#' inspired by rpart package.
#'
#' @param ... Extra arguments that would be transferred to [graphics::text()]
#' @inheritParams plot.MonoClust
#'
#' @return Labels on tree.
#'
#' @keywords internal
text_tree <- function(x,
                      which = 4,
                      digits = getOption("digits") - 2,
                      stats = TRUE,
                      abbrev,
                      cols = NULL,
                      cols.type = c("l", "p", "b"),
                      rel.loc.x = TRUE,
                      show.pval = TRUE,
                      uniform = FALSE,
                      minbranch = 0.3,
                      ...) {

  frame <- x$frame
  # These are constants that used to be arguments.
  tadj <- 0.65

  cxy <- graphics::par("cxy")
  if (!is.null(srt <- list(...)$srt) && srt == 90)
    cxy <- rev(cxy)
  xy <- plot_prep_node(x, uniform = uniform, minbranch = minbranch)

  node <- frame$number

  left_child <- match(2 * node, node)
  right_child <- match(node * 2 + 1, node)
  labels_output <- create_labels(x, abbrev = abbrev, digits = digits)
  rows <- labels_output$labels

  if (rel.loc.x) xy$x <- frame$loc + (abs(min(frame$loc))) + 1

  leaves <- frame$var == "<leaf>"
  splits <- !leaves

  left_labs <- rows[left_child[!is.na(left_child)]]
  right_labs <- rows[right_child[!is.na(right_child)]]

  # p-value display
  if (show.pval && !is.null(frame[["p.value"]])) {
    mid_labs <- frame$p.value[!is.na(frame$p.value)]
  } else {
    mid_labs <- ""
  }

  if (which == 1)
    graphics::text(xy$x[splits],
                   xy$y[splits] + tadj * cxy[2],
                   frame$var[splits], ...)
  else {
    if (which == 2 | which == 4) {
      graphics::text(xy$x[splits],
                     xy$y[splits] + tadj * cxy[2],
                     left_labs,
                     pos = 2, ...)
      # Show p-value
      if (!is.null(frame[["p.value"]]))
        graphics::text(xy$x[splits],
                       xy$y[splits] - tadj * cxy[2],
                       paste("p =", mid_labs), ...)
    }
    if (which == 3 | which == 4) {
      graphics::text(xy$x[splits],
                     xy$y[splits] + tadj * cxy[2],
                     right_labs,
                     pos = 4, ...)
      # Show p-value
      if (!is.null(frame[["p.value"]]))
        graphics::text(xy$x[splits],
                       xy$y[splits] - tadj * cxy[2],
                       paste("p =", mid_labs), ...)
    }
  }


  if (stats) {
    stat <- stringr::str_c("\n  n = ", frame$n[leaves],
                           "\n  M = ", frame$medoid[leaves])

    graphics::text(xy$x[leaves], xy$y[leaves] - tadj * cxy[2],
                   stat, adj = 0.5, ...)
  }

  # Add color bar at the bottom of the leaves
  if (!is.null(cols)) {
    if (cols.type %in% c("l", "b")) {
      graphics::rect(xy$x[leaves] - 0.05,
                     xy$y[leaves] - tadj * cxy[2] - 0.1,
                     xy$x[leaves] + 0.05,
                     xy$y[leaves] - tadj * cxy[2] - 0.08,
                     col = cols, border = NA)
    }

    if (cols.type %in% c("p", "b")) {
      graphics::points(xy$x[leaves], xy$y[leaves], pch = 16, cex = 3 *
                         graphics::par("cex"), col = cols)
    }
  }

  # Add a legend for shortened and abbreviated variable names
  if (abbrev %in% c("short", "abbreviate")) {
    varnames <- labels_output$varnames
    names <- names(varnames)
    graphics::legend(mean(xy$x), 0.9,
                     paste(varnames, names, sep   = " = "), bty = "n")
  }
}

#' Create Labels for Split Variables
#'
#' This function prints variable's labels for a `MonoClust` tree.
#'
#' @inheritParams print.MonoClust
#'
#' @return A list containing two elements:
#'   * `varnames`: A named vector of labels corresponding to variable's names
#'   (at vector names).
#'   * `labels`: Vector of labels of splitting rules to be displayed.
#' @seealso [abbreviate()]
#' @keywords internal
create_labels <- function(x, abbrev, digits = getOption("digits"), ...) {

  frame <- x$frame

  # Rename variable as indicated in abbrev
  vars <- frame$var
  uvars <- unique(vars)
  names <- uvars[uvars != "<leaf>"]

  if (abbrev == "short") {
    varnames <- stringr::str_c("V", seq_len(length(names)))
  } else if (abbrev == "abbreviate") {
    varnames <- purrr::map_chr(names, abbreviate, ...)
  } else {
    varnames <- names
  }

  names(varnames) <- names

  # Create split labels
  split_index <- which(frame$var != "<leaf>")
  lsplit <- rsplit <- character(length(split_index))

  label <- varnames[frame$var[split_index]]
  level <- frame$cut[split_index]

  # In case there is no cut information, don't show less or greater signs
  # For generalize and reuse function purposes
  if (all(is.na(level))) {
    lsplit <- rsplit <- label
  } else {
    lsplit <- paste(label, "<", round(level, digits), sep = " ")
    rsplit <- paste(label, ">=", round(level, digits), sep = " ")
  }

  node <- frame$number
  parent <- match(node %/% 2, node[split_index])
  odd <- as.logical(node %% 2)

  labels <- character(nrow(frame))
  labels[odd] <- rsplit[parent[odd]]
  labels[!odd] <- lsplit[parent[!odd]]
  labels[1] <- "root"

  return(list(varnames = varnames, labels = labels))
}
