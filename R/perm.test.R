#' Permutation Test on Monothetic Tree
#'
#' Testing the significance of each monothetic clustering split by permutation
#' methods. The original method (method 1) shuffles the observations between two
#' groups without the splitting variable. The new methods shuffle the values in
#' the splitting variable to create a new data set, then it either splits again
#' on that variable (method 2) or use all variables as the splitting candidates
#' (method 3).
#'
#' @param object The `MonoClust` object as the result of the clustering.
#' @param data The data set which is being clustered.
#' @param auto.pick Whether the algorithm stops when p-value becomes larger than
#'   `sig.val` or keeps testing and let the researcher pick the final splitting
#'   tree. Default value is `FALSE`.
#' @param sig.val Significance value to decide when to stop splitting. This
#'   option is ignored if `auto.pick = FALSE`, and is 0.05 by default when
#'   `auto.pick = TRUE`.
#' @param method Can be chosen between 1 (default), 2, or 3. See description
#'   above for the details.
#' @param rep Number of permutations required to calculate test statistic.
#' @param stat Applied to methods 2 and 3. Statistic to use, choosing between
#'   `"F"` (default) or `"AW"`.
#' @param bon.adj Whether to adjust for multiple testing problem using
#'  Bonferroni correction.
#' @inheritParams vegan::adonis
#'
#' @return The same `MonoClust` object with an extra column (p-value), as well
#' as the `numofclusters` object if `auto.pick = TRUE`.
#'
#' @details
#' ## Permutation Methods
#' ### Method 1: Shuffle the observations between two proposed clusters
#' The pseudo-F's calculated from the shuffles create the reference distribution
#' to find the p-value. Because the splitting variable that was chosen is
#' already the best in terms of reduction of inertia, that variable is withheld
#' from the distance matrix used in the permutation test.
#'
#' ### Method 2: Shuffle while keeping other variables fixed - ASW
#' This method shuffles the values of the splitting variables while keeping
#' other variables fixed to create a new data set, then the average silhouette
#' width (Kaufman and Rousseeuw, 1990) is used as the measure of separation
#' between the two new clusters and is calculated to create the reference
#' distribution.
#'
#' ### Method 3: Shuffle while keeping other variables fixed - Pseudo-F
#' Similar to the previous method but pseudo-F (as in Method 1) is used as the
#' test statistic instead of the average silhouette width.
#'
#' ## Bonferroni Correction
#' A hypothesis test occurred lower in the monothetic clustering tree could have
#' its p-value corrected for multiple tests happened before it in order to reach
#' that node. The formula is
#' \deqn{adj.p = unadj.p \times depth,}
#' with \eqn{depth} is 1 at the root node.
#'
#' @references
#' * Kaufman, L. and Rousseeuw, P. J. (1990). Finding Groups in Data: An
#' Introduction to Cluster Analysis. 1st ed. Wiley-Interscience, p. 368. isbn:
#' 978-0471735786.
#' @export
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#'
#' ruspini6sol <- MonoClust(ruspini, nclusters = 6)
#' ruspini6.p_value <- perm.test(ruspini6sol, data = ruspini, method = 1,
#'                               rep = 1000)
#' ruspini6.p_value
perm.test <- function(object, data, auto.pick = FALSE, sig.val = 0.05,
                      method = 1,
                      rep = 10000,
                      stat = c("F", "AW"),
                      bon.adj = TRUE,
                      parallel = getOption("mc.cores")) {

  if (!inherits(object, "MonoClust"))
    stop("Not a legitimate \"MonoClust\" object")
  stat <- match.arg(stat)

  frame <- object$frame
  # node <- frame$number

  # It would be better to create a jump table for reference of tree walking
  jump_table <- tibble::add_column(make_jump_table(frame),
                                   split.order = frame$split.order)
  # jump_table$right <- jump_table$left <- NA
  jump_table$p_value <- NA

  # for (i in 2:nrow(jump_table)) {
  #   parent.node <- jump_table$node[i]%/%2
  #   is.left <- jump_table$node[i]%%2 == 0
  #
  #   if (is.left) {
  #     jump_table$left[which(jump_table$node == parent.node)] <- i
  #   } else jump_table$right[which(jump_table$node == parent.node)] <- i
  # }

  jump_table$members <- NA
  jump_table$members[1] <- paste(seq_len(nrow(data)), collapse = ",")

  # Now tracing the tree to find the cluster Trace along the split order,
  # quit when reaching the
  # end or if auto.pick, stop when p value > sig.val
  last_split <- max(jump_table$split.order[!is.na(jump_table$split.order)]) + 1
  for (i in 1:max(jump_table$split.order[!is.na(jump_table$split.order)])) {
    current <- which(jump_table$split.order == i)

    # Node
    split_var <- jump_table$bipartvar[current]
    split_value <- jump_table$cut[current]

    members <- as.numeric(strsplit(jump_table$members[current], ",")[[1]])

    members_l <- members[which(data[members, split_var] < split_value)]
    jump_table$members[jump_table$left[current]] <-
      paste(members_l, collapse = ",")

    members_r <- setdiff(members, members_l)
    jump_table$members[jump_table$right[current]] <-
      paste(members_r, collapse = ",")

    p_value_unadj <- test_split(current, members, members_l, members_r,
                                auto.pick, method, data, split_var,
                                jump_table$number[current], rep, stat,
                                parallel = parallel)

    # If Bonferroni correction is applied
    p_value <- ifelse(bon.adj, p_value_unadj * i, p_value_unadj)

    jump_table$p_value[current] <- ifelse(p_value > 1, 1,
                                          ceiling(p_value * rep) / rep)

    # If auto.pick is applied
    last_pick <- NULL
    if (auto.pick) {
      if (p_value > sig.val) {
        last_split <- i
        break
      } else if (i ==
                 max(jump_table$split.order[!is.na(jump_table$split.order)])) {
        last_split <- i + 1
      }
    }
  }


  # # Tracing the tree by Node-Left-Right algorithm assign('.Jump_Table', jump_table, envir =
  # .GlobalEnv) # assign('.Data', data, envir = .GlobalEnv) recursive.walk(1, 1:nrow(data),
  # auto.pick, sig.val, method, data = data) jump_table <- .Jump_Table rm(list =
  # c('.Jump_Table'), envir = globalenv())

  frame$p.value <- jump_table$p_value

  object$frame <- frame
  object$numofclusters <- last_split
  return(object)
}

#' Hypothesis Test at Split
#'
#' @param current Row index of the current considered split.
#' @param members,members_l,members_r Vector of the index of observations that
#'   are members of the parent node, left child node, and right child node,
#'   respectively.
#' @param split_var Splitting variable at `current` split.
#' @param node The current node index
#' @inheritParams perm.test
#'
#' @return To be filled
#' @keywords internal
test_split <- function(current, members, members_l, members_r, auto.pick,
                       method, data, split_var, node, rep, stat, parallel) {
  # Membership is consecutive because the distance matrix will be moved around
  # with members_l and
  # members_r put next to each other.
  fmem2 <- c(rep(1, length(members_l)), rep(2, length(members_r)))

  if (method == 1) {
    # Method 1: shuffling the membership, not used splitting variable
    data_temp <- data
    data_temp[, split_var] <- NULL
    distmat_reduced <- as.matrix(cluster::daisy(data.frame(data_temp)))

    # A distance matrix with observations left and right put in order
    distmat_twogroup <- distmat_reduced[c(members_l, members_r),
                                        c(members_l, members_r)]

    result <- vegan::adonis(distmat_twogroup ~ fmem2,
                            permutations = rep,
                            parallel = parallel)

    # p_value.adj <- (node %/% 2 + 1) * result$aov.tab[1,6]
    p_value <- result$aov.tab$`Pr(>F)`[1]
  } else if (method == 2 | method == 3) {
    # Method 2, 3: shuffling the splitting variable, split again on that
    # variable (method 2) or split on all variables (method 3)
    currentdata <- data[c(members_l, members_r), ]

    # Find the observed statistic
    distmat_twogroup <- as.matrix(cluster::daisy(data.frame(currentdata)))
    if (stat == "F") {
      stat_obs <- F.stat(distmat_twogroup ~ fmem2, parallel = parallel)
    } else {
      stat_obs <- fpc::cluster.stats(distmat_twogroup, fmem2)$avg.silwidth
    }

    # Find the referenced distribution Create permuted matrix
    perm <- permute::shuffleSet(nrow(currentdata),
                                control = permute::how(nperm = rep))
    # The number of permutations may be smaller than rep
    permutations <- nrow(perm)

    stat_rep <- numeric(permutations)
    # Shuffling

    for (k in 1:permutations) {
      print(k)
      currentdata[, split_var] <- currentdata[perm[k, ], split_var]

      if (method == 2) {
        cluster_rep <- MonoClust(toclust = currentdata, nclusters = 2,
                                 variables = split_var)
      } else if (method == 3) {
        cluster_rep <- MonoClust(toclust = currentdata, nclusters = 2)
      }

      distmat_rep <- cluster_rep$Dist
      fmem2_rep <- cluster_rep$Membership

      if (stat == "F") {
        # If no split is made because of minbucket, cluster membership will have
        # 1 in it. In that case,
        # F-stat = 0
        stat_rep[k] <- ifelse(1 %in% fmem2_rep, 0,
                              F.stat(distmat_rep ~ fmem2_rep))
      } else {
        stat_rep[k] <- ifelse(1 %in% fmem2_rep, 0,
                              fpc::cluster.stats(distmat_rep,
                                                 fmem2_rep)$avg.silwidth)
      }
    }

    # p_value.adj <- (node %/% 2 + 1) * sum(f.stat_rep.c >= f.stat_obs) / rep
    p_value <- (sum(stat_rep >= stat_obs) + 1 -
                 sqrt(.Machine$double.eps)) / (permutations + 1)
  }

  return(p_value)
}


# ####################################################
# Copy directly from vegan::adonis function
# Modified to stop at the F.stat
# ####################################################
#



#' Title
#'
#' @inheritParams vegan::adonis
#'
#' @return Fill in later
#'
#' @importFrom stats model.frame model.matrix
F.stat <- function(formula, data = NULL, permutations = 999, method = "bray",
                   strata = NULL, contr.unordered = "contr.sum",
                   contr.ordered = "contr.poly",
                   parallel = getOption("mc.cores"), ...) {
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
  } else if ((is.matrix(lhs) || is.data.frame(lhs)) && isSymmetric(unname(as.matrix(lhs)))) {
    dmat <- as.matrix(lhs^2)
    lhs <- as.dist(lhs)
  } else {
    dist.lhs <- as.matrix(vegan::vegdist(lhs, method = method, ...))
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
