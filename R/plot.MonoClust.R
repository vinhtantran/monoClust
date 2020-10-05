#' Plot MonoClust Splitting Rule Tree
#'
#' Print the MonoClust tree in the form of dendrogram.
#'
#' @param uniform if TRUE, uniform vertical spacing of the nodes is used; this
#'   may be less cluttered when fitting a large plot onto a page. The default is
#'   to use a non-uniform spacing proportional to the error in the fit.
#' @param margin Margins to the border of the graphics machine
#' @param which Fill in later
#' @param text Whether to print the labels on the tree.
#' @param cols Whether to use color on the labels or not. Only works when `text`
#'   is `TRUE`.
#' @param rel.loc.x Whether to use the relative distance between clusters as x
#'   coordinate of the leaves. Default is TRUE.
#' @param ... Arguments to be passed to [graphics::plot.default()] and
#'   [graphics::lines()]
#'
#' @return A plot of splitting rule.
#' @export
#'
#' @importFrom graphics axis
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#' ruspini4sol <- MonoClust(ruspini, nclusters = 4)
#' plot(ruspini4sol)
#' @inheritParams print.MonoClust
plot.MonoClust <- function(x, uniform = FALSE, margin, which, abbrev = 4,
                           text = TRUE, cols = NULL, rel.loc.x = TRUE, ...) {
  ## This function sets some defaults and changes things a bit, but is mostly
  ## a wrapper for our slightly modified version of rpart's plot function (see
  ## plots.R).

  if (!inherits(x, "MonoClust"))
    stop("Not an MonoClust object")

  if (missing(margin)) {
    margin <- c(0.12, 0.02, 0, 0.05)
  }
  if (missing(which)) {
    which <- 4
  }

  plot_tree(x, uniform = uniform, margin = margin, rel.loc.x = rel.loc.x, ...)

  ## REMOVE: Tan, 3/1/15, Remove Inertia line lines(x=c(.88,.88),y=c(0,1))

  # for(i in seq(0,1,.1)) { lines(x=c(.86,.88),y=c(i,i)) text(.73,i,i) }

  if (text) {
    text_tree(x, which = which, abbrev = abbrev, cols = cols,
                   rel.loc.x = rel.loc.x)
  }

  if (!is.null(x$circularroot$var)) {
    text(x = 1, y = 1, "Circ root")
    for (i in seq_len(length(x$circularroot$var))) {
      text(x = 1, y = 1 - i * 0.05, paste(x$terms[x$circularroot$var[i]], ": ",
                                          x$circularroot$cut[i]))
    }
  }
}

Nclustplot <- function(x, main, type, ylab, xlab, ...) {

  if (missing(main)) {
    main <- "Marginal Cluster Analysis"
  }
  if (missing(type)) {
    type <- "b"
  }
  if (missing(ylab)) {
    ylab <- "Proportion of Deviance Explained"
  }
  if (missing(xlab)) {
    xlab <- "Number of Clusters"
  }

  inds <- seq(from = 2, to = nrow(x$frame), by = 2)
  plot(inds, round((1 - as.numeric(x$frame$yval[inds])/1), digits = 2), type = type, xaxt = "n",
       ylab = ylab, xlab = xlab, main = main)
  axis(1, at = inds, labels = as.character(2 + 0:(length(inds) - 1)))

}

#' Print Labels on MonoClust Tree
#'
#' @param ... Others.
#'
#' @return
#'
#' @importFrom graphics legend
#'
#' @keywords internal
#' @inheritParams plot.MonoClust
text_tree <- function(x, abbrev, which, rel.loc.x, ...) {
  ## Set up some defaults and abbreviate, then use text.rpart.
  if (missing(which))
    which <- 3
  # REMOVE: Tan, 3/1/15, remove intertia line text(.62,.5,'Inertia Explained', srt=90)
  text_MonoClust(x, which = which, abbrev = abbrev, rel.loc.x = rel.loc.x, ...)

  if (abbrev == "L") {
    vars <- x$frame$var
    uvars <- unique(vars)
    names <- uvars[uvars != "<leaf>"]
    nums <- paste("V", 1:length(names), sep = "")
    legend(mean(max(x$frame$loc), min(x$frame$loc)), 0.9, paste(nums, names, sep = " = "), bty = "n")
  }
}

#' Title
#'
#' @param x
#' @param splits
#' @param which
#' @param label
#' @param FUN
#' @param all.leaves
#' @param pretty
#' @param digits
#' @param tadj
#' @param stats
#' @param use.n
#' @param bars
#' @param legend
#' @param xadj
#' @param yadj
#' @param bord
#' @param big.pts
#' @param abbrev
#' @param cols
#' @param rel.loc.x
#' @param ...
#'
#' @return
#'
#' @importFrom graphics legend par points rect
#' @importFrom stats na.omit
#'
#' @keywords internal
#'
#' @examples
text_MonoClust <- function (x, splits = TRUE, which = 4, label = "var", FUN = text,
                     all.leaves = FALSE, pretty = NULL, digits = getOption("digits") -
                       2, tadj = 0.65, stats = TRUE, use.n = TRUE, bars = TRUE,
                     legend = FALSE, xadj = 1, yadj = 1, bord = FALSE, big.pts = FALSE, abbrev=4,
                     cols = NULL, rel.loc.x = TRUE,...)
{
  if (!inherits(x, "rpart"))
    stop("Not legitimate rpart")
  if (!is.null(x$frame$splits))
    x <- rpconvert(x)
  frame <- x$frame
  col <- names(frame)
  method <- x$method
  ylevels <- attr(x, "ylevels")
  if (!is.null(ylevels <- attr(x, "ylevels")))
    col <- c(col, ylevels)
  if (is.na(match(label, col)))
    stop("Label must be a column label of the frame component of the tree")
  cxy <- par("cxy")
  if (!is.null(srt <- list(...)$srt) && srt == 90)
    cxy <- rev(cxy)
  xy <- rpartco(x)
  #print("here")
  #print(xy)
  node <- as.numeric(row.names(x$frame))
  is.left <- (node%%2 == 0)
  node.left <- node[is.left]
  parent <- match(node.left/2, node)
  bars <- bars & is.matrix(frame$yval2)
  text.adj <- ifelse(bars, yadj * diff(range(xy$y))/12, 0)
  if (splits) {
    left.child <- match(2 * node, node)
    right.child <- match(node * 2 + 1, node)
    rows <- labels(x, abbrev=abbrev, pretty = pretty)
    #print("start")
    #print(left.child)
    #print(right.child)
    #print(rows)
    #rows <- rows[order(x$frame$loc)]

    ## Thisworks, but aren't getting them all
    ## WORK NEEDED HERE.
    #left.child <- node[left.child]
    #right.child <- node[right.child]
    #print("asdf")
    #print(rows)
    #print(node)
    #print(left.child)
    #print(right.child)

    if (rel.loc.x) xy$x <- x$frame$loc + (abs(min(x$frame$loc))) +1
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
    #print(left.child )
    #print(rows[left.child])
    #print("one")
    #print(orders)
    #print("ord")
    #print(rows[left.child[orders]])
    #print("twp")

    #print(rows[right.child[orders]])
    #print("Get labels in right order")

    leaves <- if (all.leaves)

      rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"


    splits<-!leaves
    sorders <- x$frame$loc
    #print(sorders)

    ## REMOVE: Tan, 12/14, useless command
    # labs <- rows[-1]

    ## MODIFY: Tan, 12/14, use correct left and right nodes info
    left_labs <- rows[na.omit(left.child)]
    right_labs <- rows[na.omit(right.child)]
    # left_labs <- labs[seq(from=1, to=length(labs), by=2)]
    # right_labs <- labs[seq(from=2,to=length(labs), by=2)]

    ## ADD: Tan. 3/1/15, add p-value display
    if (!is.null(x$frame$p.value)) {
      mid_labs <- na.omit(x$frame$p.value)
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
      FUN(xy$x[splits], xy$y[splits] + tadj * cxy[2], left_labs, ...)
    else {
      if (which == 2 | which == 4) {
        FUN(xy$x[splits], xy$y[splits] + tadj * cxy[2], left_labs,
            pos = 2, ...)
        # ADD: Tan, 3/1/15, Add p-value show up
        if (!is.null(frame$p.value)) FUN(xy$x[splits], xy$y[splits] - tadj * cxy[2], paste("p =", mid_labs), ...)
      }
      if (which == 3 | which == 4)
        FUN(xy$x[splits], xy$y[splits] + tadj * cxy[2], right_labs,
            pos = 4, ...)
    }
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
    if (is.null(frame$yval2))
      stat <- x$functions$text(yval = frame$yval[leaves],
                               dev = frame$dev[leaves], wt = frame$wt[leaves],
                               ylevel = ylevels, digits = digits, n = frame$n[leaves],
                               use.n = use.n, names = frame$name[leaves], meds=frame$medoid[leaves])
    else stat <- x$functions$text(yval = frame$yval2[leaves,
    ], dev = frame$dev[leaves], wt = frame$wt[leaves],
    ylevel = ylevels, digits = digits, n = frame$n[leaves],
    use.n = use.n)
    #print(xy$x[leaves])
    #print(xy$x[leaves])

    FUN(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj,
        stat,adj = 0.5, ...)
    # ADD: Tan -- 4/23/2018. Add color bar at the bottom of the leaves
    if (!is.null(cols)) {
      rect(xy$x[leaves] - 0.05, xy$y[leaves] - tadj * cxy[2] - text.adj - 0.1,
           xy$x[leaves] + 0.15, xy$y[leaves] - tadj * cxy[2] - text.adj - 0.08,
           col = cols, border = NA)
    }
  }


  #print(stat)

  #print(node[leaves])

  ## This for even
  #print(order(-1*node[leaves]))
  #print(order(node[leaves]))


  if (bars) {
    bar.vals <- x$functions$bar(yval2 = frame$yval2)
    sub.barplot(xy$x, xy$y, bar.vals, leaves, xadj = xadj,
                yadj = yadj, bord = bord, line = TRUE, col = c("lightblue",
                                                               "blue", "darkblue"))
    rx <- range(xy$x)
    ry <- range(xy$y)
    if (!is.null(ylevels))
      bar.labs <- ylevels
    else bar.labs <- dimnames(x$y)[[2]]
    if (legend & !is.null(bar.labs))
      legend(min(xy$x) - 0.1 * rx, max(xy$y) + 0.05 * ry,
             bar.labs, col = c("lightblue", "blue", "darkblue"),
             pch = 15, bty = "n", ...)
  }
  if (big.pts)
    points(xy$x[leaves], xy$y[leaves], pch = 16, cex = 3 *
             par()$cex, col = 2:(sum(leaves) + 1))
  invisible()
  if(abbrev==0){
    for(bbb in 1:length(x$qualordered)){
      if(length(x$Catnames) >0){
        print(gsub("*~*",".",x$Catnames[bbb], fixed=TRUE))
        print(cbind(1:length(x$qualordered[[bbb]]),x$qualordered[[bbb]]))
      }
    }
  }
}
