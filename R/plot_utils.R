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
                      minbranch = 0.3, rel_loc_x = TRUE, ...) {
  # REMOVE: Tan, 10/4/2020, unnecessary args
  # if (compress & missing(nspace))
  #   nspace <- branch
  # if (!compress)
  #   nspace <- -1

  # Length of the root node
  bar <- 0.03

  dev <- grDevices::dev.cur()
  if (dev == 1)
    dev <- 2

  # REMOVE, Tan, 10/4/2020, no longer use globalenv
  # assign(paste(".rpart.parms", dev, sep = "."), list(uniform = uniform,
  #                                                    branch = branch,
  #                                                    # nspace = nspace,
  #                                                    minbranch = minbranch),
  #        envir = .GlobalEnv)
  temp <- plot_prep_node(x, uniform = uniform,
                  # nspace = nspace,
                  minbranch = minbranch)
  xx <- temp$x
  yy <- temp$y

  if (rel_loc_x)
    xx <- x$frame$loc + (abs(min(x$frame$loc))) + 1
  #print("temp")
  #print(xx)

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
  # REMOVE, Tan, 10/4/2020, no use
  # method <- tree$method
  node <- as.numeric(frame$number)
  depth <- tree_depth(node)
  is_leaf <- frame$var == "<leaf>"
  if (uniform)
    y <- (1 + max(depth) - depth) / max(depth, 4)
  else {
    y <- dev <- frame$inertia
    temp <- split(seq(node), depth)
    ## change this because we are capturing number but not row

    # Index of parent of nodes in node
    parent <- match(floor(node / 2), node)

    ## REMOVE: Tan, 12/14. This line wrongly changes the labels.
    # parent[1:3] <- c(NA,1,1)

    #parrow <- which(node == parent)
    ##changed
    sibling <- match(ifelse(node %% 2, node - 1, node + 1), node)

    ## REMOVE: Tan, 12/14. This line wrongly changes the labels.
    # sibling[1:3] <- c(NA,3,2)

    #sibrow <- which(node == sibling)

    #print(parent)
    #print(parrow)
    #print(sibling)
    #print(parent)
    #print(sibrow)
    #break;
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

  ##mine
  #left_child <- node[left_child]
  #right_child <- node[right_child]
  temp <- split(seq(node)[!is_leaf], depth[!is_leaf])
  for (i in rev(temp)) x[i] <- 0.5 * (x[left_child[i]] + x[right_child[i]])
  # REMOVE, Tan, 10/4/2020 remove compress feature
  # if (nspace < 0)
  return(list(x = x, y = y))
  # REMOVE, Tan, 10/4/2020 remove compress feature
  # compress <- function(me, depth) {
  #   lson <- me + 1
  #   x <- x
  #   if (is_leaf[lson])
  #     left <- list(left = x[lson], right = x[lson], depth = depth +
  #                    1, sons = lson)
  #   else left <- compress(me + 1, depth + 1)
  #   rson <- me + 1 + length(left$sons)
  #   if (is_leaf[rson])
  #     right <- list(left = x[rson], right = x[rson], depth = depth +
  #                     1, sons = rson)
  #   else right <- compress(rson, depth + 1)
  #   maxd <- max(left$depth, right$depth) - depth
  #   mind <- min(left$depth, right$depth) - depth
  #   slide <- min(right$left[1:mind] - left$right[1:mind]) -
  #     1
  #   if (slide > 0) {
  #     x[right$sons] <- x[right$sons] - slide
  #     x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
  #     x <<- x
  #   }
  #   else slide <- 0
  #   if (left$depth > right$depth) {
  #     templ <- left$left
  #     tempr <- left$right
  #     tempr[1:mind] <- pmax(tempr[1:mind], right$right -
  #                             slide)
  #   }
  #   else {
  #     templ <- right$left - slide
  #     tempr <- right$right - slide
  #     templ[1:mind] <- pmin(templ[1:mind], left$left)
  #   }
  #   list(left = c(x[me] - nspace * (x[me] - x[lson]), templ),
  #        right = c(x[me] - nspace * (x[me] - x[rson]), tempr),
  #        depth = maxd + depth, sons = c(me, left$sons, right$sons))
  # }
  # compress(1, 1)
  # list(x = x, y = y)
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

# clusplot.MonoClust <- function(x, col.p = as.numeric(factor(x$Membership)),
#                                lines = 0,
#                                main ="Plot of Monothetic Cluster Solution",
#                                ...) {
#   cluster::clusplot(x$Dist, clus = x$Membership, col.p = col.p,
#                     lines = lines, main = main, ...)
# }

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
                      cols_type = c("l", "p", "b"),
                      rel_loc_x = TRUE,
                      show_pval = TRUE,
                      uniform = FALSE,
                      minbranch = 0.3,
                      ...) {
  # if (!inherits(x, "MonoClust"))
  #   stop("Not a legitimate MonoClust object")
  frame <- x$frame
  # col <- colnames(frame)
  # REMOVE, 10/4/2020, no use
  # method <- x$method
  # ylevels <- attr(x, "ylevels")
  # if (!is.null(ylevels <- attr(x, "ylevels")))
  #   col <- c(col, ylevels)

  # ADD, Tan, 10/10/2020, these are constants that used to be arguments.
  tadj <- 0.65
  # label <- "var"
  # REMOVE, Tan, 10/10/2020, no need to check these
  # if (is.na(match(label, col)))
  #   stop("Label must be a column label of the frame component of the tree")
  cxy <- graphics::par("cxy")
  if (!is.null(srt <- list(...)$srt) && srt == 90)
    cxy <- rev(cxy)
  xy <- plot_prep_node(x, uniform = uniform, minbranch = minbranch)
  #print("here")
  #print(xy)
  node <- frame$number
  # node_left <- node[(node %% 2) == 0]
  # parent <- match(node_left / 2, node)
  # bars <- bars & is.matrix(frame$yval2)
  # text.adj <- ifelse(bars, yadj * diff(range(xy$y))/12, 0)
  text.adj <- 0
  left_child <- match(2 * node, node)
  right_child <- match(node * 2 + 1, node)
  labels_output <- create_labels(x, abbrev = abbrev, digits = digits)
  rows <- labels_output$labels
  #print("start")
  #print(left_child)
  #print(right_child)
  #print(rows)
  #rows <- rows[order(x$frame$loc)]

  ## Thisworks, but aren't getting them all
  ## WORK NEEDED HERE.
  #left_child <- node[left_child]
  #right_child <- node[right_child]
  #print("asdf")
  #print(rows)
  #print(node)
  #print(left_child)
  #print(right_child)

  if (rel_loc_x) xy$x <- frame$loc + (abs(min(frame$loc))) + 1
  ## Find name for split



  #print(xy$x)
  #print(xy$y)
  #print("Find Split")
  #print(splits)
  #print(xy$x)
  #print(xy$y)
  #orders <- order(x$frame$loc)
  #print(orders)
  #print(xy$x[orders])
  #print(rows)

  #rows <- rows[-1]
  #print("lchild")
  #print(left_child )
  #print(rows[left_child])
  #print("one")
  #print(orders)
  #print("ord")
  #print(rows[left_child[orders]])
  #print("twp")

  #print(rows[right_child[orders]])
  #print("Get labels in right order")

  # leaves <- if (all.leaves)
  #
  #   rep(TRUE, nrow(frame))
  # else frame$var == "<leaf>"
  leaves <- frame$var == "<leaf>"


  splits <- !leaves
  # sorders <- x$frame$loc
  #print(sorders)

  ## REMOVE: Tan, 12/14, useless command
  # labs <- rows[-1]

  ## MODIFY: Tan, 12/14, use correct left and right nodes info
  left_labs <- rows[left_child[!is.na(left_child)]]
  right_labs <- rows[right_child[!is.na(right_child)]]
  # left_labs <- labs[seq(from=1, to=length(labs), by=2)]
  # right_labs <- labs[seq(from=2,to=length(labs), by=2)]

  ## ADD: Tan. 3/1/15, add p-value display
  if (show_pval && !is.null(frame[["p.value"]])) {
    mid_labs <- frame$p.value[!is.na(frame$p.value)]
  } else {
    mid_labs <- ""
  }


  #print(rows[-1])
  #print(splits)
  #print(locs <- splits*sorders)
  #print(locs2 <- locs[locs !=0 ])

  #print(left_labs[order(locs2)])
  #print(right_labs[order(locs2)])


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
      # ADD: Tan, 3/1/15, Add p-value show up
      if (!is.null(frame[["p.value"]]))
        graphics::text(xy$x[splits],
                       xy$y[splits] - tadj * cxy[2],
                       paste("p =", mid_labs), ...)
    }
    if (which == 3 | which == 4)
      graphics::text(xy$x[splits],
                     xy$y[splits] + tadj * cxy[2],
                     right_labs,
                     pos = 4, ...)
  }



  #print(node[leaves])
  #print("Right HERE")
  #print(leaves)
  #print(leaves[order(x$frame$loc)])





  #print(1:sum(splits)*2 - 1)

  ## 6 splits
  #print(splits)
  #print((sorders*splits)[1:sum(splits)])
  #print(rows[-1])


  #leaves <- rev(leaves)
  if (stats) {
    # if (is.null(frame$yval2))
    stat <- stringr::str_c("\n  n = ", frame$n[leaves],
                           "\n  M = ", frame$medoid[leaves])
    #
    # x$functions$text(yval = frame$yval[leaves],
    #                          dev = frame$dev[leaves], wt = frame$wt[leaves],
    #                          # TODO
    #                          ylevel = ylevels,
    #                          digits = digits, n = frame$n[leaves],
    #                          use.n = use.n, names = frame$name[leaves],
    #                          meds=frame$medoid[leaves])
    # else stat <- x$functions$text(yval = frame$yval2[leaves,
    # ], dev = frame$dev[leaves], wt = frame$wt[leaves],
    # # TODO
    # ylevel = ylevels,
    # digits = digits, n = frame$n[leaves],
    # use.n = use.n)
    #print(xy$x[leaves])
    #print(xy$x[leaves])

    graphics::text(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj,
        stat, adj = 0.5, ...)
  }

  # ADD: Tan -- 4/23/2018. Add color bar at the bottom of the leaves
  if (!is.null(cols)) {
    if (cols_type %in% c("l", "b")) {
      graphics::rect(xy$x[leaves] - 0.05,
                     xy$y[leaves] - tadj * cxy[2] - text.adj - 0.1,
                     xy$x[leaves] + 0.05,
                     xy$y[leaves] - tadj * cxy[2] - text.adj - 0.08,
                     col = cols, border = NA)
    }

    if (cols_type %in% c("p", "b")) {
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


  #print(stat)

  #print(node[leaves])

  ## This for even
  #print(order(-1*node[leaves]))
  #print(order(node[leaves]))


  # if (bars) {
  #   bar.vals <- x$functions$bar(yval2 = frame$yval2)
  #   sub.barplot(xy$x, xy$y, bar.vals, leaves, xadj = xadj,
  #               yadj = yadj, bord = bord, line = TRUE,
  #               col = c("lightblue",
  #               "blue", "darkblue"))
  #   rx <- range(xy$x)
  #   ry <- range(xy$y)
  #   # if (!is.null(ylevels))
  #   #   bar.labs <- ylevels
  #   # else
  #   bar.labs <- dimnames(x$y)[[2]]
  #   if (legend & !is.null(bar.labs))
  #     graphics::legend(min(xy$x) - 0.1 * rx, max(xy$y) + 0.05 * ry,
  #            bar.labs, col = c("lightblue", "blue", "darkblue"),
  #            pch = 15, bty = "n", ...)
  # }

  return(invisible(x))
  # if(abbrev==0) {
  #   for(bbb in 1:length(x$qualordered)) {
  #     if(length(x$Catnames) >0) {
  #       print(gsub("*~*",".",x$Catnames[bbb], fixed=TRUE))
  #       print(cbind(1:length(x$qualordered[[bbb]]),x$qualordered[[bbb]]))
  #     }
  #   }
  # }
}

# Nclustplot <- function(x, main, type, ylab, xlab, ...) {
#
#   if (missing(main)) {
#     main <- "Marginal Cluster Analysis"
#   }
#   if (missing(type)) {
#     type <- "b"
#   }
#   if (missing(ylab)) {
#     ylab <- "Proportion of Deviance Explained"
#   }
#   if (missing(xlab)) {
#     xlab <- "Number of Clusters"
#   }
#
#   inds <- seq(from = 2, to = nrow(x$frame), by = 2)
#   plot(inds, round((1 - as.numeric(x$frame$yval[inds]) / 1), digits = 2),
#        type = type, xaxt = "n",
#        ylab = ylab, xlab = xlab, main = main)
#   graphics::axis(1, at = inds, labels = as.character(2 + 0:(length(inds) - 1)))
#
# }
