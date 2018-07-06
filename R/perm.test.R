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
#' @param rep Number of permutations to calculate test statistic.
#'
#' @return The same MonoClust object with an updated frame with one extra
#' column (p-value), and the numofclusters chosen if auto.pick is TRUE
#' @importFrom cluster daisy
#' @export
#'
#' @examples EMPTY
perm.test <- function(object, data, auto.pick = FALSE, sig.val = 0.05,
                      method = 1, rep = 1000){

  # if(getRversion() >= "2.15.1")  utils::globalVariables(c(".Jump_Table",
  #                                                        ".Data"),
  #                                                      add = FALSE)

  if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")

  frame <- object$frame
  node <- as.numeric(row.names(frame))

  # It would be better to create a jump table for reference of tree walking
  jump.table <- cbind(node, frame[,c("bipartvar", "cut", "split.order")])
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

  jump.table$members <- NA
  jump.table$members[1] <- paste(1:nrow(data), collapse=",")

  # Now tracing the tree to find the cluster
  # Trace along the split order, quit when reaching the end or if auto.pick,
  # stop when p value > sig.val
  last.split <- max(na.omit(jump.table$split.order)) + 1
  for (i in 1:max(na.omit(jump.table$split.order))) {
    current <- which(jump.table$split.order == i)

    # Node
    split.var <- jump.table$bipartvar[current]
    split.value <- jump.table$cut[current]

    members <- as.numeric(strsplit(jump.table$members[current], ",")[[1]])

    members.L <- members[which(data[members, split.var] < split.value)]
    jump.table$members[jump.table$left[current]] <-
      paste(members.L, collapse=",")

    members.R <- setdiff(members, members.L)
    jump.table$members[jump.table$right[current]] <-
      paste(members.R, collapse=",")

    p.value.unadj <- test.split(current, members, members.L, members.R, auto.pick,
                          method, data, split.var, jump.table$node[current],
                          rep)
    p.value <- p.value.unadj * i

    jump.table$p.value[current] <- ifelse(p.value > 1, 1, p.value)

    if (auto.pick && (p.value > sig.val)) {
      last.split <- i
      break
    }
  }


  # # Tracing the tree by Node-Left-Right algorithm
  # assign(".Jump_Table", jump.table, envir = .GlobalEnv)
  # # assign(".Data", data, envir = .GlobalEnv)
  # recursive.walk(1, 1:nrow(data), auto.pick, sig.val, method,
  #                fulldata = data)
  #
  # jump.table <- .Jump_Table
  # rm(list = c(".Jump_Table"), envir = globalenv())

  frame$p.value <- jump.table$p.value
  if (auto.pick)
    object$numofclusters <- last.split

  object$frame <- frame
  return(object)
}

test.split <- function(current, members, members.L, members.R,
                       auto.pick, method, fulldata, split.var,
                       node, rep) {
  REP <- rep

  # Membership is consecutive because the distance matrix will be moved around
  # with members.L and members.R put next to each other.
  fmem2 <- factor(c(rep(1, length(members.L)), rep(2, length(members.R))))

  if (method == 1) {
    # Method 1: shuffling the membership
    data.temp <- fulldata
    data.temp[,split.var] <- NULL
    distmat.reduced <- as.matrix(daisy(data.frame(data.temp)))

    # A distance matrix with observations left and right put in order
    dist.mat.twogroup <- distmat.reduced[c(members.L, members.R),c(members.L, members.R)]

    result <- adonis(dist.mat.twogroup ~ fmem2)

    # pvalue.adj <- (node %/% 2 + 1) * result$aov.tab[1,6]
    pvalue.adj <- result$aov.tab[1,6]
  } else if (method == 2) {
    # Method 2: shuffling the splitting variable, split again on that variable
    currentdata <- fulldata[c(members.L, members.R), ]
    dist.mat.twogroup <- as.matrix(daisy(data.frame(currentdata)))

    f.stat.obs <- F.stat(dist.mat.twogroup ~ fmem2)

    f.stat.rep.c <- numeric(REP)
    # Shuffling
    for (k in 1:REP) {
      shuffled.splitting.var <- sample(currentdata[,split.var],
                                       size = nrow(currentdata),
                                       replace = FALSE)
      currentdata[,split.var] <- shuffled.splitting.var

      cluster.rep.constrained <- MonoClust(toclust = currentdata, nclusters = 2,
                                           variables = split.var)
      dist.mat.rep.c <- cluster.rep.constrained$Dist
      fmem2.rep.c <- cluster.rep.constrained$Membership

      # If no split is made because of minbucket, cluster membership will have
      # 1 in it. In that case, F-stat = 0
      f.stat.rep.c[k] <- ifelse(1 %in% fmem2.rep.c,
                                0,
                                F.stat(dist.mat.rep.c ~ fmem2.rep.c))
    }

    # pvalue.adj <- (node %/% 2 + 1) * sum(f.stat.rep.c >= f.stat.obs) / REP
    pvalue.adj <- sum(f.stat.rep.c >= f.stat.obs) / REP
  } else if (method == 3) {
    # Method 3: shuffling the splitting variable, clustering on all variables
    currentdata <- fulldata[c(members.L, members.R), ]
    dist.mat.twogroup <- as.matrix(daisy(data.frame(currentdata)))

    f.stat.obs <- F.stat(dist.mat.twogroup ~ fmem2)

    f.stat.rep.u <- numeric(REP)
    # Shuffling
    for (k in 1:REP) {
      shuffled.splitting.var <- sample(currentdata[,split.var],
                                       size = nrow(currentdata),
                                       replace = FALSE)
      currentdata[,split.var] <- shuffled.splitting.var

      cluster.rep.unconstrained <- MonoClust(toclust = currentdata, nclusters = 2)
      dist.mat.rep.u <- cluster.rep.unconstrained$Dist
      fmem2.rep.u <- cluster.rep.unconstrained$Membership
      # If no split is made because of minbucket, cluster membership will have
      # 1 in it. In that case, F-stat = 0
      f.stat.rep.u[k] <- ifelse(1 %in% fmem2.rep.u,
                                0,
                                F.stat(dist.mat.rep.u ~ fmem2.rep.u))
    }

    # pvalue.adj <- (node %/% 2 + 1) * sum(f.stat.rep.u >= f.stat.obs) / REP
    pvalue.adj <- sum(f.stat.rep.u >= f.stat.obs) / REP
  }

  return(pvalue.adj)
}


# recursive.walk <- function(current, members, auto.pick, sig.val, method,
#                            fulldata) {
#   REP <- 100
#
#   # Check stopping at leaf node
#   if (is.na(.Jump_Table$bipartvar[current])) return(0)
#
#   # Node
#   split.var <- .Jump_Table$bipartvar[current]
#   split.value <- .Jump_Table$cut[current]
#
#   # Take observations in that branch out and split into Left and Right
#   members.L <- members[which(fulldata[members, split.var] < split.value)]
#   members.R <- setdiff(members, members.L)
#
#   # Membership is consecutive because the distance matrix will be moved around
#   # with members.L and members.R put next to each other.
#   fmem2 <- factor(c(rep(1, length(members.L)), rep(2, length(members.R))))
#
#   if (method == 1) {
#     # Method 1: shuffling the membership
#     data.temp <- fulldata
#     data.temp[,split.var] <- NULL
#     distmat.reduced <- as.matrix(daisy(data.frame(data.temp)))
#
#     # A distance matrix with observations left and right put in order
#     dist.mat.twogroup <- distmat.reduced[c(members.L, members.R),c(members.L, members.R)]
#
#     result <- adonis(dist.mat.twogroup ~ fmem2)
#     pvalue.adj <- (as.numeric(row.names(.Jump_Table[current,])) %/% 2 + 1) * result$aov.tab[1,6]
#   } else if (method == 2) {
#     # Method 2: shuffling the splitting variable, split again on that variable
#     currentdata <- fulldata[c(members.L, members.R), ]
#     dist.mat.twogroup <- as.matrix(daisy(data.frame(currentdata)))
#
#     f.stat.obs <- F.stat(dist.mat.twogroup ~ fmem2)
#
#     f.stat.rep.c <- numeric(REP)
#     # Shuffling
#     for (k in 1:REP) {
#       shuffled.splitting.var <- sample(currentdata[,split.var],
#                                        size = nrow(currentdata),
#                                        replace = FALSE)
#       currentdata[,split.var] <- shuffled.splitting.var
#
#       cluster.rep.constrained <- MonoClust(toclust = currentdata, nclusters = 2,
#                                            variables = split.var)
#       dist.mat.rep.c <- cluster.rep.constrained$Dist
#       fmem2.rep.c <- cluster.rep.constrained$Membership
#       f.stat.rep.c[k] <- F.stat(dist.mat.rep.c ~ fmem2.rep.c)
#     }
#
#     pvalue.adj <- sum(f.stat.rep.c < f.stat.obs) / REP
#   } else if (method == 3) {
#     # Method 3: shuffling the splitting variable, clustering on all variables
#     dist.mat.twogroup <- as.matrix(daisy(data.frame(data.temp)))
#
#     f.stat.obs <- F.stat(dist.mat.twogroup ~ fmem2)
#
#     f.stat.rep.u <- numeric(REP)
#     # Shuffling
#     for (k in 1:REP) {
#       shuffled.splitting.var <- sample(data.temp[,split.var],
#                                        size = nrow(data.temp),
#                                        replace = FALSE)
#       data.temp[,split.var] <- shuffled.splitting.var
#
#       cluster.rep.unconstrained <- MonoClust(toclust = data.temp, nclusters = 2)
#       dist.mat.rep.u <- cluster.rep.unconstrained$Dist
#       fmem2.rep.u <- cluster.rep.unconstrained$Membership
#       f.stat.rep.u[k] <- F.stat(dist.mat.rep.u ~ fmem2.rep.u)
#     }
#
#     pvalue.adj <- sum(f.stat.rep.u < f.stat.obs) / REP
#   }
#
#   .Jump_Table$p.value[current] <<- ifelse(pvalue.adj > 1, 1, pvalue.adj)
#
#   if (auto.pick && (pvalue.adj > sig.val)) return(0)
#   # Left
#   recursive.walk(.Jump_Table$left[current], members.L, auto.pick, sig.val, method, fulldata)
#
#   # Right
#   recursive.walk(.Jump_Table$right[current], members.R, auto.pick, sig.val, method, fulldata)
# }

#' ####################################################
#' Copy directly from adonis function
#' Modified to stop at the F.stat
#' ####################################################
F.stat <- function(formula, data = NULL, permutations = 999, method = "bray",
                   strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly",
                   parallel = getOption("mc.cores"), ...)
{
  TOL <- 1e-07
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  qrhs <- qr(rhs)
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1)
    stop("right-hand-side of formula has no usable terms")
  H.s <- lapply(2:length(u.grps), function(j) {
    Xj <- rhs[, grps %in% u.grps[1:j]]
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    tcrossprod(Q[, 1:qrX$rank])
  })
  if (inherits(lhs, "dist")) {
    if (any(lhs < -TOL))
      stop("dissimilarities must be non-negative")
    dmat <- as.matrix(lhs^2)
  }
  else if ((is.matrix(lhs) || is.data.frame(lhs)) && isSymmetric(unname(as.matrix(lhs)))) {
    dmat <- as.matrix(lhs^2)
    lhs <- as.dist(lhs)
  }
  else {
    dist.lhs <- as.matrix(vegdist(lhs, method = method,
                                  ...))
    dmat <- dist.lhs^2
  }
  n <- nrow(dmat)
  G <- -sweep(dmat, 1, rowMeans(dmat))/2
  SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
  SS.Exp.each <- c(SS.Exp.comb - c(0, SS.Exp.comb[-nterms]))
  H.snterm <- H.s[[nterms]]
  tIH.snterm <- t(diag(n) - H.snterm)

  SS.Res <- sum(G * tIH.snterm)
  df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
  df.Res <- n - qrhs$rank

  F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
  return(F.Mod)
}
#' END COPY ########################################################
