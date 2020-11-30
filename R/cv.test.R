#' Cross-Validation Test on MonoClust
#'
#' Perform cross-validation test for different different number of clusters of
#' Monothetic Clustering.
#'
#' @param data Data set to be partitioned.
#' @param fold Number of folds (k). `fold = 1` is the special case, when the
#' function performs a Leave-One-Out Cross-Validation (LOOCV).
#' @param minnodes Minimum number of clusters to be checked.
#' @param maxnodes Maximum number of clusters to be checked.
#' @param ncores Number of CPU cores on the current host. When set to NULL,
#'   all available cores are used.
#' @param ... Other parameters transferred to [MonoClust()].
#'
#' @details
#' The \eqn{k}-fold cross-validation randomly partitions data into \eqn{k}
#' subsets with equal (or close to equal) sizes. \eqn{k - 1} subsets are used as
#' the training data set to create a tree with a desired number of leaves and
#' the other subset is used as validation data set to evaluate the predictive
#' performance of the trained tree. The process repeats for each subset as the
#' validating set (\eqn{m = 1, \ldots, k}) and the mean squared difference,
#' \deqn{MSE_m=\frac{1}{n_m} \sum_{q=1}^Q\sum_{i \in m} d^2_{euc}(y_{iq},
#' \hat{y}_{(-i)q}),}
#' is calculated, where \eqn{\hat{y}_{(-i)q}} is the cluster mean on the
#' variable
#' \eqn{q} of the cluster created by the training data where the observed value,
#' \eqn{y_{iq}}, of the validation data set will fall into, and
#' \eqn{d^2_{euc}(y_{iq}, \hat{y}_{(-i)q})} is the squared Euclidean distance
#' (dissimilarity) between two observations at variable $q$. This process is
#' repeated for the $k$ subsets of the data set and the average of these test
#' errors is the cross-validation-based estimate of the mean squared error of
#' predicting a new observation,
#' \deqn{CV_K = \overline{MSE} = \frac{1}{M} \sum_{m=1}^M MSE_m.}
#'
#' @note This function supports parallel processing with [foreach::foreach()].
#'   It distributes MonoClust calls to processes.
#'
#' @return A `MonoClust.cv` class containing a data frame of mean sum of square
#'   error and its standard deviation.
#' @seealso [plot.cv.MonoClust()], [MonoClust()], [predict.MonoClust()]
#' @export
#'
#' @examples
#' \donttest{
#' library(cluster)
#' data(ruspini)
#'
#' # Leave-one-out cross-validation
#' cv.test(ruspini, fold = 1, minnodes = 2, maxnodes = 4)
#'
#' # 5-fold cross-validation
#' cv.test(ruspini, fold = 5, minnodes = 2, maxnodes = 4)
#' }
cv.test <- function(data, fold = 10L, minnodes = 2L, maxnodes = 10L,
                    ncores = 1L, ...) {

  if (!is.data.frame(data)) {
    stop("\"data\" must be a data frame.")
  }

  if (!is.null(ncores)){
    if (!is.numeric(ncores)) {
      stop("\"ncores\" should be either NULL or a positive integer.")
    }
    if (ncores < 1) {
      stop("\"ncores\" should be > 1.")
    }
  }

  if (is.null(ncores))
    ncores <- parallel::detectCores() - 1

  if (ncores > 1) {
    # Initiate processes
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    `%dodo%` <- foreach::`%dopar%`
  } else {
    `%dodo%` <- foreach::`%do%`
  }


  num_obs <- nrow(data)
  sse_t <- vector("list", maxnodes - minnodes + 1)
  # LOOCV
  if (fold == 1) {
    for (k in minnodes:maxnodes) {
      sse_i <-
        foreach::foreach(iter = seq_len(num_obs),
                         .combine = "c",
                         .inorder = FALSE,
                         .packages = c("monoClust")) %dodo% {
                           out <- MonoClust(data[-iter, ], nclusters = k, ...)
                           pred <- predict.MonoClust(out,
                                                     newdata = data[iter, ],
                                                     type = "centroid")
                           return(sum((data[iter, ] - pred[, -1])^2))
                         }
      sse_t[[k - minnodes + 1]] <- c(ncluster = k,
                                     MSE = mean(sse_i),
                                     `Std. Dev.` = stats::sd(sse_i))
    }

    ret <- list(cv = dplyr::bind_rows(sse_t),
                cv.type = "Leave-one-out Cross-validation")
  } else {
    # sse_t is a list of vector, later will be bound into a tibble
    sse_t <- vector("list", maxnodes - minnodes + 1)
    index <- rep(1:fold, num_obs %/% fold + 1)
    random_list <- sample(index, num_obs, replace = FALSE)

    for (k in minnodes:maxnodes) {
      sse_i <-
        foreach::foreach(iter = 1:fold,
                         .combine = "c",
                         .inorder = FALSE,
                         .packages = c("monoClust")) %dodo% {
                           train_set <- data[-which(random_list == iter), ]
                           test_set <- data[which(random_list == iter), ]
                           train_tree <- MonoClust(train_set, nclusters = k)
                           pred <- predict.MonoClust(train_tree,
                                                     test_set,
                                                     type = "centroid")
                           return(sum((test_set - pred[, -1])^2))
                         }
      sse_t[[k - minnodes + 1]] <- c(ncluster = k,
                                     MSE = mean(sse_i),
                                     `Std. Dev.` = stats::sd(sse_i))
    }
    ret <- list(cv = dplyr::bind_rows(sse_t),
                cv.type = paste0(fold, "-fold Cross-validation"))
  }

  if (ncores > 1) {
    # Stop processes
    parallel::stopCluster(cl)
  }

  class(ret) <- "cv.MonoClust"
  return(ret)
}
