cv.test <- function(data, fold = 10, minnodes = 2, maxnodes = 10, ...) {

    SSET=numeric(0)
    num.obs <- dim(data)[1]
    # LOOCV
    if (fold == 1)
    {
        for (k in minnodes:maxnodes){
            SSEi=numeric(0)
            fullc<- MonoClust(da1, nclusters = k)
            # pRsq=c(pRsq,adonis(dist(da1)~factor(fullc$Membership),permutations=1)$aov.tab[1,5])
            for (i in seq(length(da1[,1]))){
                out <- MonoClust(da1[-i,], nclusters = k)
                predict(out,newdata=da1[i,])
                SSEi=c(SSEi,sum((da1[i,]-predict(out,newdata=da1[i,])[,-1])^2))
            }

            SSET=c(SSET,sum(SSEi))

        }
    }
    else if (fold > 1) {
        index <- rep(1:fold, num.obs %/% fold + 1)
        random.list <- sample(index, num.obs, replace=FALSE)

        for (k in minnodes:maxnodes) {
            SSEi=numeric(0)
            for (i in 1:fold) {
                train.set <- data[-which(random.list == i),]
                test.set <- data[which(random.list == i),]
                train.tree <- MonoClust(train.set, nclusters = k)
                SSEi=c(SSEi, sum((test.set-predict(train.tree, test.set)[,-1])^2))
            }
            SSET=rbind(SSET,c(mean(SSEi), sd(SSEi)))
            colnames(SSET) <- c("MSE", "Std. Dev.")
        }
    }

    SSET

}
