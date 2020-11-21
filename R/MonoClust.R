#' Monothetic Clustering
#'
#' Creates a MonoClust object after partitioning the data set using Monothetic
#' Clustering.
#'
#' @param toclust Data set as a data frame.
#' @param cir.var Index or name of the circular variable in the data set.
#' @param variables List of variables selected for clustering procedure. It
#'   could be a vector of variable indexes, or a vector of variable names.
#' @param distmethod Distance method to use with the data set. The default value
#'   is the Euclidean distance but Gower will be used if there is circular
#'   variable (`cir.var` is specified). Transfer to [cluster::daisy()].
#' @param digits Significant decimal number printed in the output.
#' @param nclusters Number of clusters created. Default is 2.
#' @param minsplit The minimum number of observations that must exist in a node
#'   in order for a split to be attempted. Default is 5.
#' @param minbucket The minimum number of observations in any terminal leaf
#'   node. Default is `minsplit/3`.
#' @param ncores Number of CPU cores on the current host. If greater than 1,
#'   parallel processing with [foreach::foreach()] is used to distribute cut
#'   search on variables to processes. When set to NULL, all available cores are
#'   used.
#'
#' @return A `MonoClust` object. See [`MonoClust.object`].
#' @importFrom rlang .data
#' @export
#'
#' @references
#' 1. Chavent, M. (1998). A monothetic clustering method. Pattern Recognition
#' Letters, 19(11), 989-996. \doi{10.1016/S0167-8655(98)00087-7}.
#' 2. Tran, T. V. (2019). Monothetic Cluster Analysis with Extensions to
#' Circular and Functional Data. Montana State University - Bozeman.
#'
#' @examples
#' # Very simple data set
#' library(cluster)
#' data(ruspini)
#' ruspini4sol <- MonoClust(ruspini, nclusters = 4)
#' ruspini4sol
#'
#' # data with circular variable
#' library(monoClust)
#' data(wind_sensit_2007)
#'
#' # Use a small data set
#' set.seed(12345)
#' wind_reduced <- wind_sensit_2007[sample.int(nrow(wind_sensit_2007), 10), ]
#' circular_wind <- MonoClust(wind_reduced, cir.var = 3, nclusters = 2)
#' circular_wind
MonoClust <- function(toclust,
                      cir.var = NULL,
                      variables = NULL,
                      distmethod = c("euclidean", "manhattan", "gower"),
                      digits = getOption("digits"),
                      nclusters = 2L,
                      minsplit = 5L,
                      minbucket = round(minsplit / 3),
                      ncores = 1L) {

  if (!is.data.frame(toclust)) {
    stop("\"toclust\" must be a data frame.")
  }

  toclust <- dplyr::as_tibble(toclust)

  if (minbucket >= minsplit) {
    stop("\"minbucket\" must be less than \"minsplit\".")
  }

  # Check that no categorical variables existed
  if (!all(purrr::map_lgl(toclust, is.numeric))) {
    stop("Function does not support categorical variables yet.")
  }

  # Argument checking (variables)
  if (!is.null(variables)) {
    if (!is.vector(variables)) {
      stop("Variables need to be a vector of variable names or indices.")
    }
    toclust1 <- toclust[variables]
    if (!is.data.frame(toclust1)) {
      stop("Undefined variables selected.")
    }
    if (variables %in% colnames(toclust)) {
      variables <- which(colnames(toclust) == variables)
    }
  } else {
    variables <- seq_len(ncol(toclust))
  }

  # Check the distmethod
  if (is.null(distmethod)) {
    # If there is a circular variable in the data set, use Gower's
    # distance unless otherwise specified.
    distmethod <- ifelse(!is.null(cir.var), "gower", "euclidean")
  } else {
    distmethod <- match.arg(distmethod)
  }

  # Argument checking (circular variables)
  if (!is.null(cir.var)) {
    if (!is.vector(cir.var)) {
      stop("Circular variables need to be a vector of variable names or
           indices.")
    } else if (length(cir.var) > 1) {
      stop("MonoClust supports only one circular variable in the data set.")
    }
    if (!is.data.frame(toclust[cir.var])) {
      stop("Undefined variables selected.")
    }
  }

  # Impute missing values if mice is installed. Emit warning.
  if (any(is.na(toclust))) {
    if (requireNamespace("mice", quietly = TRUE)) {
      imputed <- mice::mice(toclust)
      message(strwrap(c("Data contain missing values. mice() is used for
                        imputation. See mice() help page for more details.
                        Missing cells per column:", imputed$nmis),
                      prefix = "\n", initial = ""))
      toclust <- mice::complete(imputed)
    }
    else
      stop(strwrap("Data contain NA. Install \"mice\" package and rerun this
                   function to automatically impute the missing value(s).",
                   prefix = "\n", initial = ""))
  }

  if (!is.null(ncores)){
    if (!is.numeric(ncores)) {
      stop("\"ncores\" should be either NULL or a positive integer")
    }
    if (ncores < 1) {
      stop("\"ncores\" should be > 1")
    }
  }

  if (is.null(ncores))
    ncores <- parallel::detectCores() - 1

  bestcircsplit <- NULL
  if (!is.null(cir.var)) {

    next_value <- purrr::map_dfc(toclust[cir.var], find_closest)
    bestcircsplit <- list(hour = 0, minute = 0, inertia = 0)
    distmat_circ <- as.matrix(circ_dist(toclust[, cir.var]) * length(cir.var))
    inertia_circ <- inertia_calc(distmat_circ)

    variable <- toclust[[cir.var]]
    min_value <- min(variable)
    max_value <- max(variable)
    min_inertia <- Inf
    # Check all starting value to find the best split
    # The best split will be recorded in output with the pivot is hour
    while (min_value < max_value) {

      # Set the hour hand
      variable_shift <- variable %cd-% min_value

      new_toclust <- dplyr::tibble(variable_shift)
      out <-
        checkem(data = new_toclust,
                cuts = purrr::map_dfc(new_toclust, find_closest),
                frame = new_node(number            = 1L,
                                 var               = "<leaf>",
                                 n                 = nrow(new_toclust),
                                 inertia           = inertia_circ,
                                 inertia_explained = 0,
                                 medoid            = 1,
                                 loc               = 0.1,
                                 split.order       = 0L),
                cloc = rep(1, nrow(new_toclust)),
                dist = distmat_circ,
                variables = 1,
                minsplit = minsplit,
                minbucket = minbucket,
                split_order = 0,
                ncores = ncores)

      cut <- out$frame$cut[1]
      inertia <- sum(out$frame[out$frame$var == "<leaf>", "inertia"])

      if (min_inertia > inertia) {
        min_inertia <- inertia
        bestcircsplit <- list(hour = min_value,
                              minute = cut %cd+% min_value,
                              intertia = inertia)
      }
      # Increase min_value to the next higher value
      min_value <- dplyr::pull(next_value[which(variable == min_value), 1])[1]
    }

    # Shift the circular variable to the hour pivot. That turns the circular
    # to linear variable. The first best split would be the minute split,
    # which was already found, together with hour, to be the best arcs. Will
    # shift back later by modifying cluster_frame
    toclust[, cir.var] <- variable %cd-% bestcircsplit$hour

  }

  cuts <- purrr::map_dfc(toclust, find_closest)
  nobs <- nrow(toclust)

  # Add metric argument to daisy and circular distance
  if (!is.null(cir.var)) {
    # Have to split because R coerce the distance matrix to one value when
    # using ifelse with a 0.
    # distmat1: distance matrix of pure quantitative variables
    if (ncol(toclust[, -cir.var]) != 0L) {
      distmat1 <- cluster::daisy(toclust[, -cir.var], metric = distmethod) *
        ncol(toclust[, -cir.var])
    } else {
      distmat1 <- 0
    }
    # distmat2: distance matrix of circular variables
    # circ_var accepts multiple circular variables
    distmat2 <- circ_dist(toclust[, cir.var]) * ncol(toclust[, cir.var])

    # distmat0: combined from 1 and 2 as the mean distances
    distmat0 <- (distmat1 + distmat2) / ncol(toclust)
  } else
    distmat0 <- cluster::daisy(toclust, metric = distmethod)

  distmat <- as.matrix(distmat0)

  members <- seq_len(nobs)
  c_loc <- rep(1, nobs)

  cluster_frame <-
    new_node(number            = 1L,
             var               = "<leaf>",
             n                 = nobs,
             inertia           = inertia_calc(distmat[members, members]),
             inertia_explained = 0,
             medoid            = medoid(members, distmat),
             loc               = 0.1,
             split.order       = 0L)

  split_order <- 1L
  done_running <- FALSE
  # This loop runs until we have nclusters, have exhausted our observations or
  # run into our minbucket/minsplit restrictions.
  while (sum(cluster_frame$var == "<leaf>") < nclusters && !done_running) {

    checkem_ret <- checkem(toclust, cuts, cluster_frame, c_loc, distmat,
                           variables,
                           minsplit, minbucket, split_order, ncores)
    split_order <- split_order + 1L

    # Use c_loc because it is only ran in splitter
    if (!identical(c_loc, checkem_ret$cloc)) {
      cluster_frame <- checkem_ret$frame
    } else {
      done_running <- TRUE
    }

    c_loc <- checkem_ret$cloc
  }

  # Modify the Cluster_frame to shift the circular variable's cut back to the
  # original values
  if (!is.null(cir.var)) {
    cir_pos <- which(cluster_frame$var == colnames(toclust)[cir.var])
    cluster_frame$cut[cir_pos] <-
      cluster_frame$cut[cir_pos] %cd+% bestcircsplit$hour
  }

  # Add medoids of each cluster
  medoids <- cluster_frame$medoid[cluster_frame$var == "<leaf>"]
  names(medoids) <- cluster_frame$number[cluster_frame$var == "<leaf>"]

  # For display purpose, all -99 is turned to NA
  cluster_frame <- tibble::add_column(
    replace(dplyr::select(cluster_frame, -.data$alt),
            dplyr::select(cluster_frame, -.data$alt) == -99,
            NA),
    dplyr::select(cluster_frame, .data$alt)
  )

  # Make variables corresponding to leaves NA
  cluster_frame[cluster_frame$var == "<leaf>",
                c("bipartsplitrow", "bipartsplitcol")] <- NA

  # Remove all rows from alt data frame if it is a leave
  cluster_frame$alt <-
    purrr::map_if(cluster_frame$alt,
                  cluster_frame$var == "<leaf>",
                  ~ slice(.x, 0))

  # Whether there exists an alternate splitting route
  alt <- sum(purrr::map_int(frame$alt, nrow)) > 0

  MonoClust_obj <-
    list(frame = cluster_frame,
         membership = c_loc,
         dist = distmat,
         # Add terms to keep track of variables name
         terms = colnames(toclust),
         # Centroids info, for prediction
         centroids = centroid(toclust, cluster_frame, c_loc),
         medoids = medoids,
         alt = alt,
         circularroot = list(var = cir.var, cut = bestcircsplit$hour)
    )

  class(MonoClust_obj) <- "MonoClust"

  return(MonoClust_obj)
}

#' Split Function
#'
#' Given the Cluster's frame's row position to split at `split_row`, this
#' function performs the split, calculate all necessary information for the
#' splitting tree and cluster memberships.
#'
#' @param split_row The row index in frame that would be split on.
#' @inheritParams checkem
#'
#' @return Updated `frame` and `cloc` saved in a list.
#'
#' @keywords internal
splitter <- function(data, cuts, split_row, frame, cloc, dist,
                     split_order = 0L) {

  node_number <- frame$number[split_row]
  mems <- which(cloc == node_number)
  split <- c(frame$bipartsplitrow[split_row],
             frame$bipartsplitcol[split_row])

  # Extract data and cuts to split
  datamems <- data[mems, ]
  cutsmems <- cuts[mems, ]

  # Split into data rows of lower half (A) and upper half (B)
  mems_a <- mems[which(datamems[[split[2]]] < cutsmems[[split[1], split[2]]])]
  mems_b <- setdiff(mems, mems_a)

  mid_cutpoint <- mean(c(datamems[[split[1], split[2]]],
                         cutsmems[[split[1], split[2]]]))

  # Make the new clusters.
  node_number_a <- node_number * 2L
  node_number_b <- node_number * 2L + 1L

  cloc[mems_a] <- node_number_a
  cloc[mems_b] <- node_number_b

  variable_name <- colnames(data)[frame$bipartsplitcol[split_row]]

  frame$var[split_row]         <- variable_name
  frame$cut[split_row]         <- mid_cutpoint
  frame$split.order[split_row] <- split_order

  # New cluster 1
  node_a <-
    new_node(
      number  = node_number_a,
      var     = "<leaf>",
      n       = length(mems_a),
      inertia = inertia_calc(dist[mems_a, mems_a]),
      medoid  = medoid(mems_a, dist),
      loc     = frame$loc[split_row] - 1 / nrow(frame)
    )

  # New cluster 2
  node_b <-
    new_node(
      number  = node_number_b,
      var     = "<leaf>",
      n       = length(mems_b),
      inertia = inertia_calc(dist[mems_b, mems_b]),
      medoid  = medoid(mems_b, dist),
      loc     = frame$loc[split_row] + 1 / nrow(frame)
    )

  # Insert two new rows right after split row
  frame <- tibble::add_row(frame,
                           tibble::add_row(node_a, node_b),
                           .after = split_row)

  # This has to be updated last because it needs leaf nodes list
  # See Chavent (2007) for definition. Basically,
  # 1 - (sum(current inertia)/inertia[1])
  frame$inertia_explained[split_row] <-
    1 - sum(frame$inertia[frame$var == "<leaf>"]) / frame$inertia[1]

  return(list(frame = frame, cloc = cloc))
}

#' Find the Best Split
#'
#' Find the best split in terms of reduction in inertia for the transferred
#' node, indicate by row. Find the terminal node with the greatest change in
#' inertia and bi-partition it.
#'
#' @param frame_row One row of the split tree as data frame.
#' @inheritParams checkem
#'
#' @return The updated `frame_row` with the next split updated.
#' @keywords internal
find_split <- function(data, cuts, frame_row, cloc, dist, variables, minsplit,
                       minbucket, ncores) {

  node_number <- frame_row$number
  mems <- which(cloc == node_number)
  inertiap <- frame_row$inertia

  if (inertiap == 0L || frame_row$n < minsplit || frame_row$n == 1L) {
    frame_row$bipartsplitrow <- 0L
    return(frame_row)
  }

  # Subset the data and cut matricies
  datamems <- data[mems, variables]
  cutsmems <- cuts[mems, variables]

  # For each possible cut, calculate the inertia. This is where using a
  # discrete optimization algorithm would help a lot.
  mult_inertia <- function(i, datamems, cutsmems, dist, mems) {
    data_col <- dplyr::pull(datamems, i)
    cuts_col <- dplyr::pull(cutsmems, i)

    new_inertia <- purrr::map_dbl(cuts_col, function(x) {
      mems_a <- mems[which(data_col < x)]
      mems_b <- setdiff(mems, mems_a)
      ifelse(length(mems_a) * length(mems_b) == 0L,
             NA,
             inertia_calc(dist[mems_a, mems_a]) +
               inertia_calc(dist[mems_b, mems_b]))
    })

    return(new_inertia)
  }

  # Initiate processes
  `%op%` <- foreach::`%do%`

  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    `%op%` <- foreach::`%dopar%`
  }

  bycol <-
    foreach::foreach(
      i = seq_len(ncol(datamems)),
      .combine = cbind) %op%
    mult_inertia(i, datamems, cutsmems, dist, mems)

  if (ncores > 1)
    # Stop processes
    parallel::stopCluster(cl)

  # Difference between current cluster and the possible splits
  vals <- inertiap - bycol
  # Say no difference if we have NA or infinite (happens when no split is
  # possible)
  vals[!is.finite(vals) | is.na(vals)] <- 0L

  # This is the best split
  maxval <- max(vals)

  # This is the maximum inertia change index
  ind <- which((inertiap - bycol) == maxval, arr.ind = TRUE)

  # Add one more column to check minbucket
  ind_1 <- cbind(ind, minbucket = TRUE)

  for (i in seq_len(nrow(ind_1))) {
    split <- ind_1[i, ]
    left_size <- sum(datamems[, split[2L]] < cutsmems[[split[1L], split[2L]]])
    right_size <- length(mems) - left_size
    if (left_size < minbucket | right_size < minbucket)
      ind_1[i, ncol(ind_1)] <- FALSE
  }

  # Remove all row doesn't satisfy minbucket and make sure output is always a
  # matrix even though it has only one row
  ind_1 <- matrix(ind_1[ind_1[, ncol(ind_1)] != FALSE, ],
                  ncol = ncol(ind_1))

  # If multiple splits produce the same inertia change output a warning.
  if (nrow(ind_1) > 1) {
    a <- ind_1[, -3]
    colnames(a) <- c("bipartsplitrow", "bipartsplitcol")
    frame_row$alt <- list(tibble::as_tibble(a)[-1, ])
  }

  # If there is at least one row that satisfies minbucket, pick the first one
  if (nrow(ind_1) != 0L) {
    split <- ind_1[1L, ]

    mems_a <-
      mems[
        which(datamems[, split[2L]] < cutsmems[[split[1L], split[2L]]])
      ]
    mems_b <- setdiff(mems, mems_a)

    # calculate our change in inertia
    inertiadel <- inertiap -
      inertia_calc(dist[mems_a, mems_a]) -
      inertia_calc(dist[mems_b, mems_b])

    # Update frame
    frame_row$bipartsplitrow <- split[1L]
    frame_row$bipartsplitcol <- variables[split[2L]]
    frame_row$inertiadel <- inertiadel
  } else
    # Otherwise, stop as a leaf
    frame_row$bipartsplitrow <- 0L

  return(frame_row)
}

#' First Gate Function
#'
#' This function checks what are available nodes to split and then call
#' `find_split()` on each node, then decide which node creates best split, and
#' call `splitter()` to perform the split.
#'
#' @param data Original data set.
#' @param cuts Cuts data set, which has the next higher value of each variable
#'   in the original data set.
#' @param frame The split tree transferred as data frame.
#' @param cloc Vector of current cluster membership.
#' @param dist Distance matrix of all observations in the data.
#' exported function yet. Vector of 1 for all observations.
#' @param minsplit The minimum number of observations that must exist in a node
#'   in order for a split to be attempted.
#' @param split_order The control argument to see how many split has been done.
#' @param ncores Number of CPU cores on the current host.
#' @inheritParams MonoClust
#'
#' @return It is not supposed to return anything because global environment was
#'   used. However, if there is nothing left to split, it returns 0 to tell the
#'   caller to stop running the loop.
#' @keywords internal
checkem <- function(data, cuts, frame, cloc, dist, variables, minsplit,
                    minbucket, split_order, ncores) {

  # Current terminal nodes
  candidates <- which(frame$var == "<leaf>" &
                        frame$bipartsplitrow == -99L)
  # Split the best one. Return to nada which never gets output
  frame[candidates, ] <-
    purrr::map_dfr(candidates,
                   ~ find_split(data, cuts, frame[.x, ], cloc, dist, variables,
                                minsplit, minbucket, ncores))


  # See which ones are left
  candidates2 <- which(frame$var == "<leaf>" & frame$bipartsplitrow != 0L)

  # if there is something, run. Otherwise, frame and cloc are not updated, cloc
  # is used to check if done running
  if (length(candidates2) > 0L) {
    # Find the best inertia change of all that are possible
    split_row <- candidates2[which.max(frame$inertiadel[candidates2])]

    # Make new clusters from that cluster
    splitter_ret <- splitter(data, cuts, split_row, frame, cloc, dist,
                             split_order)

    frame <- splitter_ret$frame
    cloc <- splitter_ret$cloc
  }

  return(list(frame = frame, cloc = cloc))
}
