#' Title
#'
#' @param object
#' @param newdata
#' @param na.action
#' @param type
#' @param ...
#'
#' @return
#' @export
#'
#' @importFrom stats na.pass
#'
#' @examples
predict.MonoClust <- function(object, newdata = list(), na.action = na.pass, type = "mean", ...) {

    ####### DELETE: For debugging purpose ########### require(cluster) data(ruspini) source('MonoClust.R')
    ####### object <- MonoClust(ruspini, nclusters = 4) newdata <- list(x = c(48, 90), y = c(93, 70))
    ####### na.action <- na.pass

    if (!inherits(object, "rpart"))
        stop("Not a legitimate \"rpart\" object")

    where <- if (missing(newdata))
        object$Membership else {
        frame <- object$frame
        new.data <- data.frame(newdata)
        # Check to see if the newdata has the same structure as the object
        Terms <- object$terms
        new.terms <- colnames(new.data)
        if (!length(intersect(Terms, new.terms)) == length(Terms))
            stop("The new data set does not have the same variables as the original data set")

        node <- as.numeric(row.names(frame))
        depth <- tree_depth(node)

        # It would be better to create a jump table for reference of tree walking
        jump.table <- cbind(node, frame[, c("bipartvar", "cut")])
        jump.table$right <- jump.table$left <- NA

        if (nrow(jump.table) >= 2) {
            for (i in 2:nrow(jump.table)) {
                parent.node <- jump.table$node[i]%/%2
                is.left <- jump.table$node[i]%%2 == 0

                if (is.left) {
                  jump.table$left[which(jump.table$node == parent.node)] <- i
                } else jump.table$right[which(jump.table$node == parent.node)] <- i
            }
        } else jump.table$bipartvar <- jump.table$cut <- NA

        # Now tracing the tree to find the cluster
        apply(new.data, 1, tree.walk, jump.table = jump.table)
    }



    mtype <- missing(type)
    type <- match.arg(type)
    ylevels <- attr(object, "ylevels")
    nclass <- length(ylevels)
    if (mtype && nclass > 0L)
        type <- "prob"
    if (type == "mean" && is.null(frame$yval2)) {
        centroids <- data.frame(object$centroids)
        ret <- centroids[sapply(where, function(x) which(x == centroids[, 1])), ]
    }
    # } else if (type == 'matrix') { pred <- frame$yval2[where, ] dimnames(pred) <-
    # list(names(where), NULL) } else if (type == 'class' && nclass > 0L) { if (length(ylevels) ==
    # 0L) stop('type 'class' is only appropriate for classification') pred <-
    # factor(ylevels[frame$yval[where]], levels = ylevels) names(pred) <- names(where) } else if
    # (type == 'prob' && nclass > 0L) { pred <- frame$yval2[where, 1L + nclass + 1L:nclass, drop =
    # FALSE] dimnames(pred) <- list(names(where), ylevels) } else stop('Invalid prediction for
    # \'rpart\' object')

    # Expand out the missing values in the result But only if operating on the original dataset if
    # (missing(newdata) && !is.null(object$na.action)) pred <- naresid(object$na.action, pred)
    ret

}

tree.walk <- function(jump.table, new.point) {
    current.node <- 1
    repeat {
        # Check stopping at leaf node
        if (is.na(jump.table$bipartvar[current.node]))
            break

        # If it's not a leaf node, trace
        var <- jump.table$bipartvar[current.node]
        value <- new.point[var]

        current.node <- as.numeric(ifelse(value < jump.table$cut[current.node], jump.table$left[current.node],
            jump.table$right[current.node]))
    }

    # Extract the leaf info out. TODO: Output the necessary info.  QUESTION: Do we need to
    # re-calculate any properties? Like medoid?
    node <- as.numeric(row.names(jump.table)[current.node])
    names(node) <- current.node
    node
}
