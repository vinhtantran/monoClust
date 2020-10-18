#' Permutation Test on Monothetic Tree
#'
#' Testing the significance of each monothetic clustering split by permutation
#' methods. The "simple-withhold" method (`"sw"`) shuffles the observations
#' between two groups without the splitting variable. The other two methods
#' shuffle the values in the splitting variable to create a new data set, then
#' it either splits again on that variable ("resplit-limit", `"rl"`) or use all
#' variables as the splitting candidates ("resplit-nolimit", `"rn"`).
#'
#' @param object The `MonoClust` object as the result of the clustering.
#' @param data The data set which is being clustered.
#' @param auto.pick Whether the algorithm stops when p-value becomes larger than
#'   `sig.val` or keeps testing and let the researcher pick the final splitting
#'   tree. Default value is `FALSE`.
#' @param sig.val Significance value to decide when to stop splitting. This
#'   option is ignored if `auto.pick = FALSE`, and is 0.05 by default when
#'   `auto.pick = TRUE`.
#' @param method Can be chosen between `sw` (simple-withhold, default), `rl`
#'   (resplit-limit), or `rn` (resplit-nolimit). See Details.
#' @param rep Number of permutations required to calculate test statistic.
#' @param stat Statistic to use. Choosing between `"f"` (Calinski-Harabaz's
#'   pseudo-F (Calinski and Harabasz, 1974)) or `"aw"` (Average silhoutte width
#'   by Rousseeuw (1987)).
#' @param bon.adj Whether to adjust for multiple testing problem using
#'   Bonferroni correction.
#'
#' @return The same `MonoClust` object with an extra column (p-value), as well
#' as the `numofclusters` object if `auto.pick = TRUE`.
#'
#' @details
#' ## Permutation Methods
#' ### Simple-Withhold: Shuffle the observations between two proposed clusters
#' The `stat` calculated from the shuffles create the reference distribution
#' to find the p-value. Because the splitting variable that was chosen is
#' already the best in terms of reduction of inertia, that variable is withheld
#' from the distance matrix used in the permutation test.
#'
#' ### Resplit-Limit: Shuffle splitting variable, split again on that variable
#' This method shuffles the values of the splitting variables while keeping
#' other variables fixed to create a new data set, then the chosen `stat` is
#' calculated for each rep to compare with the observed `stat`.
#'
#' ### Resplit-Nolimit: Shuffle splitting variable, split on all variables
#' Similar to Method 2 but all variables are splitting candidates.
#'
#' ## Bonferroni Correction
#' A hypothesis test occurred lower in the monothetic clustering tree could have
#' its p-value corrected for multiple tests happened before it in order to reach
#' that node. The formula is
#' \deqn{adj.p = unadj.p \times depth,}
#' with \eqn{depth} is 1 at the root node.
#'
#' @note This function uses [foreach::foreach()] to facilitate parallel
#'   processing. It distributes reps to processes.
#'
#' @references
#' Calinski, T. and Harabasz, J (1974). "A dendrite method for cluster
#' analysis". en. In: *Communications in Statistics* 3.1, pp. 1-27.
#'
#' Rousseeuw, P. J. (1987). "Silhouettes: A graphical aid to the interpretation
#' and validation of cluster analysis". In: *Journal of Computational and
#' Applied Mathematics* 20, pp. 53-65. ISSN: 03770427. DOI:
#' 10.1016/0377-0427(87) 90125-7.
#'
#' @export
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#'
#' \donttest{
#' ruspini6sol <- MonoClust(ruspini, nclusters = 6)
#' ruspini6.p_value <- perm.test(ruspini6sol, data = ruspini, method = "sw",
#'                               rep = 1000)
#' ruspini6.p_value
#' }
#'
#' \dontrun{
#' # Multiple processing via doParallel
#' library(doParallel)
#'
#' cl <- makePSOCKcluster(5)
#' registerDoParallel(cl)
#'
#' ruspini6.p_value <- perm.test(ruspini6sol, data = ruspini, method = "sw",
#'                               rep = 1000)
#'
#' stopCluster(cl)
#' }
perm.test <- function(object, data, auto.pick = FALSE, sig.val = 0.05,
                      method = c("sw", "rl", "rn"),
                      rep = 1000,
                      stat = c("f", "aw"),
                      bon.adj = TRUE) {

  if (!inherits(object, "MonoClust"))
    stop("Not a legitimate \"MonoClust\" object")
  method <- match.arg(method)
  stat <- match.arg(stat)

  frame <- object$frame

  # It would be better to create a jump table for reference of tree walking
  jump_table <- tibble::add_column(make_jump_table(frame),
                                   split.order = frame$split.order)
  jump_table$p_value <- NA

  jump_table$members <- NA
  jump_table$members[1] <- paste(seq_len(nrow(data)), collapse = ",")

  # Now tracing the tree to find the cluster Trace along the split order,
  # quit when reaching the end or if auto.pick, stop when p value > sig.val
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

    p_value_unadj <- test_split(members_l, members_r, method, data, split_var,
                                rep, stat)

    # If Bonferroni correction is applied
    p_value <- ifelse(bon.adj, p_value_unadj * i, p_value_unadj)

    jump_table$p_value[current] <- ifelse(p_value > 1, 1,
                                          ceiling(p_value * rep) / rep)

    # If auto.pick is applied
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

  frame$p.value <- jump_table$p_value

  object$frame <- frame
  object$numofclusters <- last_split
  return(object)
}

#' Hypothesis Test at Split
#'
#' @param members_l,members_r Vector of the index of observations that
#'   are members of the left child node and the right child node, respectively.
#' @param split_var Splitting variable at current split.
#' @inheritParams perm.test
#'
#' @return p-value of the test
#' @keywords internal
test_split <- function(members_l, members_r, method, data, split_var, rep,
                       stat) {
  # Membership is consecutive because the distance matrix will be moved around
  # with members_l and members_r put next to each other.
  fmem2 <- c(rep(1, length(members_l)), rep(2, length(members_r)))
  `%op%` <- get_oper(foreach::getDoParWorkers() > 1)

  if (method == "sw") {
    # Simple-withhold: shuffling the membership, not used splitting variable
    data_reduced <- data[colnames(data) != split_var]
  } else {
    data_reduced <- data
  }

  # A distance matrix with observations left and right put in order
  current_data <- data_reduced[c(members_l, members_r), , drop = FALSE]

  # Find the observed statistic
  distmat_twogroup <- as.matrix(cluster::daisy(current_data))
  stat_obs_l <- cluster_stats(distmat_twogroup, fmem2)

  # Create permuted matrix
  perm <- permute::shuffleSet(nrow(current_data),
                              control = permute::how(nperm = rep))
  # The number of permutations may be smaller than rep
  permutations <- nrow(perm)

  # Shuffling
  stat_rep_l <-
    if (method == "sw") {
      foreach::foreach(k = 1:permutations,
                       .inorder = FALSE) %op% {
                         distmat_rep <- distmat_twogroup[perm[k, ], perm[k, ]]
                         return(cluster_stats(distmat_rep, fmem2))
                       }
    } else if (method == "rl") {
      foreach::foreach(k = 1:permutations,
                       .inorder = FALSE,
                       .packages = c("monoClust")) %op% {
                         current_data[, split_var] <- current_data[perm[k, ],
                                                                   split_var]
                         # Resplit-limit: limit the splitting variables
                         cluster_rep <- MonoClust(toclust = current_data,
                                                  nclusters = 2,
                                                  variables = split_var)
                         distmat_rep <- cluster_rep$Dist
                         fmem2_rep <- cluster_rep$Membership
                         return(cluster_stats(distmat_rep, fmem2_rep))
                       }
    } else if (method == "rn") {
      foreach::foreach(k = 1:permutations,
                       .inorder = FALSE,
                       .packages = c("monoClust")) %op% {
                         current_data[, split_var] <- current_data[perm[k, ],
                                                                   split_var]
                         # Resplit-nolimit: splitting all variables
                         cluster_rep <- MonoClust(toclust = current_data,
                                                  nclusters = 2)
                         distmat_rep <- cluster_rep$Dist
                         fmem2_rep <- cluster_rep$Membership
                         return(cluster_stats(distmat_rep, fmem2_rep))
                       }
    }

  stat_rep_tbl <- dplyr::bind_rows(stat_rep_l)

  if (stat == "F") {
    stat_obs <- stat_obs_l$f_stat
    stat_rep <- stat_rep_tbl$f_stat
  } else {
    stat_obs <- stat_obs_l$asw
    stat_rep <- stat_rep_tbl$asw
  }

  p_value <- (sum(stat_rep >= stat_obs) + 1 -
                sqrt(.Machine$double.eps)) / (permutations + 1)

  return(p_value)
}


#' Cluster Statistics Calculation
#'
#' Calinski-Harabaz's pseudo-F (Calinski and Harabasz, 1974) and Average
#' silhoutte width (Rousseeuw, 1987) calculation.
#'
#' @param d Distance object (as generated by [dist()]) or a distance matrix
#'   between cases.
#' @param clustering Integer vector of length of the number of cases, which
#'   indicates a clustering. The clusters have to be numbered from 1 to the
#'   number of clusters.
#'
#' @return
#' \describe{
#'   \item{f_stat}{Calinski-Harabasz's pseudo-F.}
#'   \item{asw}{Average silhouette width.}
#' }
#'
#' @references
#' Caliński, T. and Harabasz, J (1974). "A dendrite method for cluster
#' analysis". en. In: *Communications in Statistics* 3.1, pp. 1–27.
#'
#' Rousseeuw, P. J. (1987). "Silhouettes: A graphical aid to the interpretation
#' and validation of cluster analysis". In: *Journal of Computational and
#' Applied Mathematics* 20, pp. 53–65. ISSN: 03770427. DOI:
#' 10.1016/0377-0427(87) 90125-7.
#'
#' @seealso [cluster::silhouette()]
#'
#' @keywords internal
cluster_stats <- function(d, clustering) {

  d <- as.dist(d)
  nclust <- max(clustering)
  n <- length(clustering)

  # Total cluster sum of squares
  total_ss <- sum(d^2) / n

  # Within cluster sum of squares
  dmat <- as.matrix(d)
  within_ss <- 0
  for (i in 1:nclust) {
    d_i <- as.dist(dmat[clustering == i, clustering == i])

    if (i <= nclust) {
      within_ss <- within_ss + sum(d_i^2) / sum(clustering == i)
    }
  }

  # Between cluster sum of squares
  between_ss <- total_ss - within_ss

  f_stat <- between_ss * (n - nclust) / (within_ss * (nclust - 1))

  # Average silhouette width
  sii <- cluster::silhouette(clustering, dmatrix = dmat)
  sc <- summary(sii)
  asw <- sc$avg.width

  return(list(f_stat = f_stat, asw = asw))
}
