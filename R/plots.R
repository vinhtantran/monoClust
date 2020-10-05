#' Plot the splitting tree. Inherited from rpart::plot_tree
#'
#' @param x MonoClust object
#' @param uniform see Documentation of rpart package
#' @param branch see Documentation of rpart package
#' @param compress see Documentation of rpart package
#' @param nspace see Documentation of rpart package
#' @param margin see Documentation of rpart package
#' @param minbranch see Documentation of rpart package
#' @param bar see Documentation of rpart package
#' @param rel.loc.x Whether to use the relative distance between clusters as x
#'   coordinate of the leaves. Default is TRUE.
#' @param ... more args
#'
#' @return Plot of tree
#'
#' @importFrom graphics lines
#'
#' @examples
#' @keywords internal
plot_tree <- function(x, uniform = FALSE, branch = 1, compress = FALSE, nspace,
                      margin = 0, minbranch = 0.3, bar = 0.03, rel.loc.x = TRUE,
                      ...)
{
  if (compress & missing(nspace))
    nspace <- branch
  if (!compress)
    nspace <- -1
  dev <- grDevices::dev.cur()
  if (dev == 1)
    dev <- 2
  assign(paste(".rpart.parms", dev, sep = "."), list(uniform = uniform,
                                                     branch = branch, nspace = nspace, minbranch = minbranch),
         envir = .GlobalEnv)
  temp <- rpartco(x, uniform = uniform, nspace = nspace,
                  minbranch = minbranch)
  xx <- temp$x
  yy <- temp$y

  if (rel.loc.x) xx <- x$frame$loc + (abs(min(x$frame$loc))) +1
  #print("temp")
  #print(xx)

  temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
  temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
  plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "",
       ...)
  node <- as.numeric(row.names(x$frame))
  temp <- plot_branch(xx, yy, node, branch)
  if (branch > 0)
    lines(c(xx[1], xx[1]), c(yy[1], yy[1] + bar * diff(range(yy))),
          ...)
  lines(c(temp$x), c(temp$y))
  invisible(list(x = xx, y = yy))
}

rpartco <- function(tree, uniform = FALSE, nspace = -1, minbranch = 0.3)
{
  frame <- tree$frame
  method <- tree$method
  node <- as.numeric(frame$number)
  depth <- tree_depth(node)
  is.leaf <- frame$var == "<leaf>"
  if (uniform)
    y <- (1 + max(depth) - depth)/max(depth, 4)
  else {
    y <- dev <- frame$inertia
    temp <- split(seq(node), depth)
    ## change this because we are capturing number but not row

    parent <- match(floor(node/2), node)

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
    fudge <- minbranch * diff(range(y))/max(depth)
    for (i in temp[-1]) {

      temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
      haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
      y[i] <- y[parent[i]] - ifelse(temp2 <= fudge & haskids,
                                    fudge, temp2)
    }
    y <- y/(max(y))
  }
  x <- double(length(node))
  x[is.leaf] <- seq(sum(is.leaf))
  left.child <- match(node * 2, node)
  right.child <- match(node * 2 + 1, node)

  ##mine
  #left.child <- node[left.child]
  #right.child <- node[right.child]
  temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
  for (i in rev(temp)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
  if (nspace < 0)
    return(list(x = x, y = y))
  compress <- function(me, depth) {
    lson <- me + 1
    x <- x
    if (is.leaf[lson])
      left <- list(left = x[lson], right = x[lson], depth = depth +
                     1, sons = lson)
    else left <- compress(me + 1, depth + 1)
    rson <- me + 1 + length(left$sons)
    if (is.leaf[rson])
      right <- list(left = x[rson], right = x[rson], depth = depth +
                      1, sons = rson)
    else right <- compress(rson, depth + 1)
    maxd <- max(left$depth, right$depth) - depth
    mind <- min(left$depth, right$depth) - depth
    slide <- min(right$left[1:mind] - left$right[1:mind]) -
      1
    if (slide > 0) {
      x[right$sons] <- x[right$sons] - slide
      x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
      x <<- x
    }
    else slide <- 0
    if (left$depth > right$depth) {
      templ <- left$left
      tempr <- left$right
      tempr[1:mind] <- pmax(tempr[1:mind], right$right -
                              slide)
    }
    else {
      templ <- right$left - slide
      tempr <- right$right - slide
      templ[1:mind] <- pmin(templ[1:mind], left$left)
    }
    list(left = c(x[me] - nspace * (x[me] - x[lson]), templ),
         right = c(x[me] - nspace * (x[me] - x[rson]), tempr),
         depth = maxd + depth, sons = c(me, left$sons, right$sons))
  }
  compress(1, 1)
  list(x = x, y = y)
}

plot_branch=function (x, y, node, branch = 0)
{
  is.left <- (node%%2 == 0)
  node.left <- node[is.left]
  parent <- match(node.left/2, node)
  sibling <- match(node.left + 1, node)
  temp <- (x[sibling] - x[is.left]) * (1 - branch)/2
  xx <- rbind(x[is.left], x[is.left] + temp, x[sibling] - temp,
              x[sibling], NA)
  yy <- rbind(y[is.left], y[parent], y[parent], y[sibling],
              NA)
  list(x = xx, y = yy)
}


clusplot.MonoClust <- function(x, col.p = as.numeric(factor(x$Membership)), lines = 0, main ="Plot of Monothetic Cluster Solution",...){
  cluster::clusplot(x$Dist,clus=x$Membership,col.p=col.p,lines=lines,main=main, ...)
}

