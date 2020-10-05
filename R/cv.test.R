#' k-fold Cross validation test on a monothetic clustering result.
#'
#' @param data data set
#' @param fold number of folds (k)
#' @param minnodes the minimum number of clusters to be checked
#' @param maxnodes the maximum number of clusters to be checked
#' @param ...
#'
#' @return A list of sum of squares of difference between the predicted and true
#' values
#' @export
#'
#' @importFrom stats sd
#'
#' @examples NULL
cv.test <- function(data, fold = 10, minnodes = 2, maxnodes = 10, ...) {

    SSET <- numeric(0)
    num.obs <- dim(data)[1]
    # LOOCV
    if (fold == 1) {
        for (k in minnodes:maxnodes) {
            SSEi <- numeric(0)
            fullc <- MonoClust(data, nclusters = k)
            # pRsq=c(pRsq,adonis(dist(data)~factor(fullc$Membership),permutations=1)$aov.tab[1,5])
            for (i in seq(length(data[, 1]))) {
                out <- MonoClust(data[-i, ], nclusters = k)
                predict.MonoClust(out, newdata = data[i, ])
                SSEi <- c(SSEi, sum((data[i, ] - predict.MonoClust(out, newdata = data[i, ])[, -1])^2))
            }

            SSET <- c(SSET, sum(SSEi))

        }
    } else if (fold > 1) {
        index <- rep(1:fold, num.obs%/%fold + 1)
        random.list <- sample(index, num.obs, replace = FALSE)

        for (k in minnodes:maxnodes) {
            SSEi <- numeric(0)
            for (i in 1:fold) {
                train.set <- data[-which(random.list == i), ]
                test.set <- data[which(random.list == i), ]
                train.tree <- MonoClust(train.set, nclusters = k)
                SSEi <- c(SSEi, sum((test.set - predict.MonoClust(train.tree, test.set)[, -1])^2))
            }
            SSET <- rbind(SSET, c(mean(SSEi), sd(SSEi)))

        }
        colnames(SSET) <- c("MSE", "Std. Dev.")
        rownames(SSET) <- seq(minnodes, maxnodes, 1)
    }

    SSET

}
