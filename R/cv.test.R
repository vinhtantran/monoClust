#' Cross-Validation Test on MonoClust
#'
#' Perform cross-validation test for different different number of clusters of
#' Monothetic Clustering.
#'
#' @param data Data set to be partitioned.
#' @param fold Number of folds (k). When `k = 1`, the function performs a
#' Leave-One-Out-Cross-Validation (LOOCV). When `k > 1`, a data will be
#' returned.
#' @param minnodes Minimum number of clusters to be checked.
#' @param maxnodes Maximum number of clusters to be checked.
#' @param ... Other parameters transferred to [MonoClust()].
#'
#' @details
#' The \eqn{k}-fold cross-validation randomly partitions data into \eqn{k} subsets with
#' equal (or close to equal) sizes. \eqn{k - 1} subsets are used as the training
#' data set to create a tree with a desired number of leaves and the other
#' subset is used as validation data set to evaluate the predictive performance
#' of the trained tree. The process repeats for each subset as the validating
#' set (\eqn{m = 1, \ldots, k}) and the mean squared difference,
#' \deqn{MSE_m=\frac{1}{n_m} \sum_{q=1}^Q\sum_{i \in m} d^2_{euc}(y_{iq}, \hat{y}_{(-i)q}),}
#' is calculated, where \eqn{\hat{y}_{(-i)q}} is the cluster mean on the variable
#' \eqn{q} of the cluster created by the training data where the observed value,
#' \eqn{y_{iq}}, of the validation data set will fall into, and
#' \eqn{d^2_{euc}(y_{iq}, \hat{y}_{(-i)q})} is the squared Euclidean distance
#' (dissimilarity) between two observations at variable $q$. This process is
#' repeated for the $k$ subsets of the data set and the average of these test
#' errors is the cross-validation-based estimate of the mean squared error of
#' predicting a new observation,
#' \deqn{CV_K = \overline{MSE} = \frac{1}{M} \sum_{m=1}^M MSE_m.}
#'
#' @return A list of sum of squares of difference between the predicted and true
#'   values.
#' @seealso [cv.plot()], [MonoClust()], [predict.MonoClust()]
#' @export
#'
#' @importFrom stats sd
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#' cv.test(ruspini, minnodes = 2, maxnodes = 4)
cv.test <- function(data, fold = 10, minnodes = 2, maxnodes = 10, ...) {

  if (!is.data.frame(data)) {
    stop("\"data\" must be a data frame.")
  }

  num_obs <- nrow(data)
  # LOOCV
  if (fold == 1) {
    # sse_t is a vector for LOOCV
    sse_t <- vector("double", maxnodes - minnodes + 1)
    for (k in minnodes:maxnodes) {
      sse_i <- vector("double", num_obs)
      # fullc <- MonoClust(data, nclusters = k)
      for (i in seq_len(num_obs)) {
        out <- MonoClust(data[-i, ], nclusters = k, ...)
        predict.MonoClust(out, newdata = data[i, ])
        sse_i[i] <-
          sum((data[i, ] -
                 predict.MonoClust(out,
                                   newdata = data[i, ])[, -1])^2)
      }

      sse_t[k] <- sum(sse_i)

    }
  } else if (fold > 1) {
    # sse_t is a list of vector, later will be bound into a tibble
    sse_t <- vector("list", maxnodes - minnodes + 1)
    index <- rep(1:fold, num_obs %/% fold + 1)
    random_list <- sample(index, num_obs, replace = FALSE)

    for (k in minnodes:maxnodes) {
      sse_i <- vector("double", fold)
      for (i in 1:fold) {
        train_set <- data[-which(random_list == i), ]
        test_set <- data[which(random_list == i), ]
        train_tree <- MonoClust(train_set, nclusters = k)
        sse_i[i] <-
          sum((test_set -
                 predict.MonoClust(train_tree,
                                   test_set)[, -1])^2)
      }
      sse_t[k] <- c(mean(sse_i), sd(sse_i))

    }
    sse_t <- purrr::flatten_dfr(sse_t)
    colnames(sse_t) <- c("MSE", "Std. Dev.")
    rownames(sse_t) <- seq(minnodes, maxnodes, 1)
  }

  sse_t

}
