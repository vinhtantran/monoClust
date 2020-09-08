#' Plot the splitting tree. Inherited from rpart::plot.rpart
#'
#' @param x MonoClust object
#' @param uniform see Documentation of rpart package
#' @param branch see Documentation of rpart package
#' @param compress see Documentation of rpart package
#' @param nspace see Documentation of rpart package
#' @param margin see Documentation of rpart package
#' @param minbranch see Documentation of rpart package
#' @param bar see Documentation of rpart package
#' @param rel.loc.x Whether to use the relative distance between clusters as x coordinate of the leaves. Default is TRUE.
#' @param ... more args
#'
#' @return Plot of tree
#' @export
#'
#' @examples Blank
plot.rpart=function (x, uniform = FALSE, branch = 1, compress = FALSE, nspace,
                     margin = 0, minbranch = 0.3, bar = 0.03, rel.loc.x = TRUE,
                     ...)
{
    if (!inherits(x, "rpart"))
        stop("Not an rpart object")
    if (!is.null(x$frame$splits))
        x <- rpconvert(x)
    if (compress & missing(nspace))
        nspace <- branch
    if (!compress)
        nspace <- -1
    dev <- dev.cur()
    if (dev == 1)
        dev <- 2
    assign(paste(".rpart.parms", dev, sep = "."), list(uniform = uniform,
                                                       branch = branch, nspace = nspace, minbranch = minbranch),
           envir = .GlobalEnv)
    temp <- rpartco(x)
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
    temp <- rpart.branch(xx, yy, node, branch)
    if (branch > 0)
        lines(c(xx[1], xx[1]), c(yy[1], yy[1] + bar * diff(range(yy))),
              ...)
    lines(c(temp$x), c(temp$y))
    invisible(list(x = xx, y = yy))
}

rpartco=function (tree, parms = paste(".rpart.parms", dev.cur(), sep = "."))
{
    frame <- tree$frame
    method <- tree$method
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    is.leaf <- (frame$var == "<leaf>")
    if (exists(parms, envir = .GlobalEnv)) {
        parms <- get(parms, envir = .GlobalEnv)
        uniform <- parms$uniform
        nspace <- parms$nspace
        minbranch <- parms$minbranch
    }
    else {
        uniform <- FALSE
        nspace <- -1
        minbranch <- 0.3
    }
    if (uniform)
        y <- (1 + max(depth) - depth)/max(depth, 4)
    else {
        y <- dev <- frame$dev
        temp <- split(seq(node), depth)
        ## change this because we are capturing number but not row

        parent <- match(floor(node/2), node)

        ## REMOVE: Tan, 12/14. This line wrongly changes the labels.
        # parent[1:3] <- c(NA,1,1)

        #parrow <- which(node == parent)
        ##changed
        sibling <- match(ifelse(node%%2, node - 1, node + 1), node)

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

tree.depth=function (nodes)
{
    depth <- floor(log(nodes, base = 2) + 1e-07)
    as.vector(depth - min(depth))
}

rpart.branch=function (x, y, node, branch)
{
    if (missing(branch)) {
        if (exists(parms <- paste(".rpart.parms", dev.cur(),
                                  sep = "."), envir = .GlobalEnv)) {
            parms <- get(parms, envir = .GlobalEnv)
            branch <- parms$branch
        }
        else branch <- 0
    }
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
text.rpart=function (x, splits = TRUE, which = 4, label = "var", FUN = text,
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


clusplot.MonoClust <- function(x, col.p = as.numeric(factor(x$Membership)), lines = 0, main ="Plot of Monothetic Cluster Solution",...){
    cluster::clusplot(x$Dist,clus=x$Membership,col.p=col.p,lines=lines,main=main, ...)
}

