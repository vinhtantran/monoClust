#' Title
#'
#' @param data
#' @param fold
#' @param minnodes
#' @param maxnodes
#' @param ...
#'
#' @return
#' @export
#'
#' @importFrom stats sd
#'
#' @examples
cv.test.mse <- function(data, fold = 10, minnodes = 2, maxnodes = 10, ...) {

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
            MSEi <- numeric(0)
            for (i in 1:fold) {
                train.set <- data[-which(random.list == i), ]
                test.set <- data[which(random.list == i), ]
                train.tree <- MonoClust(train.set, nclusters = k)
                MSEi <- c(MSEi, mean((test.set - predict.MonoClust(train.tree, test.set)[, -1])^2))
            }
            SSET <- rbind(SSET, c(mean(MSEi), sd(MSEi)))
            colnames(SSET) <- c("MSE", "Std. Dev.")
        }
    }

    SSET

}
