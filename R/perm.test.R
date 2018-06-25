#' Testing the significance of each monothetic clustering split by permutation
#' methods. The original method (method 1) shuffles the observations between two
#' groups without the splitting variable. The new methods shuffle the values in
#' the splitting variable to create a new data set, then it either splits again
#' on that variable (method 2) or use all variables as the splitting candidates
#' (method 3).
#'
#' @param object The MonoClust object as the result of the clustering.
#' @param data The data set which is being clustered.
#' @param auto.pick Whether the algorithm stops when p-value becomes larger than
#' sig.val or keeps testing and let the researcher pick the final splitting
#' tree. Default value is FALSE.
#' @param sig.val Significance value to decide when to stop splitting. This
#' option is ignored if auto.pick is FALSE.
#' @param method Can be chosen between 1 (default), 2, or 3. See description
#' above the details.
#'
#' @return The same MonoClust object with an updated frame with one extra
#' column (p-value), and the numofclusters chosen if auto.pick is TRUE
#' @importFrom cluster daisy
#' @export
#'
#' @examples EMPTY
perm.test <- function(object, data, auto.pick = FALSE, sig.val = 0.05,
                      method = 1){

  if(getRversion() >= "2.15.1")  utils::globalVariables(c(".Jump_Table",
                                                          ".Data"))

  if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")

  frame <- object$frame
  node <- as.numeric(row.names(frame))

  # It would be better to create a jump table for reference of tree walking
  jump.table <- cbind(node, frame[,c("bipartvar", "cut")])
  jump.table$right <- jump.table$left <- NA
  jump.table$p.value <- NA

  for (i in 2:nrow(jump.table)) {
    parent.node <- jump.table$node[i] %/% 2
    is.left <- jump.table$node[i] %% 2 == 0

    if (is.left) {
      jump.table$left[which(jump.table$node == parent.node)] <- i
    } else
      jump.table$right[which(jump.table$node == parent.node)] <- i
  }

  # Now tracing the tree to find the cluster
  # Tracing the tree by Node-Left-Right algorithm
  assign(".Jump_Table", jump.table, envir = .GlobalEnv)
  assign(".Data", data, envir = .GlobalEnv)
  recursive.walk(1, 1:nrow(data), auto.pick, sig.val, method)

  jump.table <- .Jump_Table
  rm(list = c(".Jump_Table", ".Data"), envir = globalenv())
  frame$p.value <- jump.table$p.value
  if (auto.pick) {
    a <- frame[which(frame$p.value > sig.val),]
    last.split <- min(a$split.order)
    object$numofclusters <- last.split
  }
  object$frame <- frame
  return(object)
}

recursive.walk <- function(current, members, auto.pick, sig.val, method) {
  REP <- 100

    # Check stopping at leaf node
    if (is.na(.Jump_Table$bipartvar[current])) return(0)

  # Node
  split.var <- .Jump_Table$bipartvar[current]
  split.value <- .Jump_Table$cut[current]

  data.temp <- .Data

  if (method == 1) {
    # Method 1: shuffling the membership
    data.temp[,split.var] <- NULL
    distmat.reduced <- as.matrix(daisy(data.frame(data.temp)))
    members.L <- members[which(.Data[members, split.var] < split.value)]
    members.R <- setdiff(members, members.L)

    dist.mat.twogroup <- distmat.reduced[c(members.L, members.R),c(members.L, members.R)]
    fmem2 <- factor(c(rep(1, length(members.L)), rep(2, length(members.R))))
    result <- adonis(dist.mat.twogroup ~ fmem2)
    pvalue.adj <- (as.numeric(row.names(.Jump_Table[current,])) %/% 2 + 1) * result$aov.tab[1,6]
  } else if (method == 2) {
    # Method 2: shuffling the splitting variable, split again on that variable
    dist.mat.twogroup <- as.matrix(daisy(data.frame(data.temp)))

    f.stat.obs <- F.stat(dist.mat.twogroup ~ fmem2)

    # Shuffling
    for (k in 1:REP) {
      shuffled.splitting.var <- sample(data.temp[,split.var],
                                       size = nrow(data.temp),
                                       replace = FALSE)
      data.temp[,split.var] <- shuffled.splitting.var

      cluster.rep.constrained <- MonoClust(toclust = data.temp, nclusters = 2,
                                           variables = split.var)
      dist.mat.rep.c <- cluster.rep.constrained$Dist
      fmem2.rep.c <- cluster.rep.constrained$Membership
      f.stat.rep.c[k] <- F.stat(dist.mat.rep.c ~ fmem2.rep.c)
    }

    pvalue.adj <- sum(f.stat.rep.c < f.stat.obs) / REP
  } else if (method == 3) {
    # Method 3: shuffling the splitting variable, clustering on all variables
    dist.mat.twogroup <- as.matrix(daisy(data.frame(data.temp)))

    f.stat.obs <- F.stat(dist.mat.twogroup ~ fmem2)

    # Shuffling
    for (k in 1:REP) {
      shuffled.splitting.var <- sample(data.temp[,split.var],
                                       size = nrow(data.temp),
                                       replace = FALSE)
      data.temp[,split.var] <- shuffled.splitting.var

      cluster.rep.unconstrained <- MonoClust(toclust = data.temp, nclusters = 2)
      dist.mat.rep.u <- cluster.rep.unconstrained$Dist
      fmem2.rep.u <- cluster.rep.unconstrained$Membership
      f.stat.rep.u[k] <- F.stat(dist.mat.rep.u ~ fmem2.rep.u)
    }

    pvalue.adj <- sum(f.stat.rep.u < f.stat.obs) / REP
  }

  .Jump_Table$p.value[current] <<- ifelse(pvalue.adj > 1, 1, pvalue.adj)

  if (auto.pick && (pvalue.adj > sig.val)) return(0)
  # Left
  recursive.walk(.Jump_Table$left[current], members.L, auto.pick, sig.val)

  # Right
  recursive.walk(.Jump_Table$right[current], members.R, auto.pick, sig.val)
}
