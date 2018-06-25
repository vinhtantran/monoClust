data(arctic)
summary(arctic)

obj.lst <- NULL
for (i in 1978:2014) {
  extent <- arctic[arctic$Year == i, "Extent"]
  count <-  nrow(arctic[arctic$Year == i,])

  date.list <- with(arctic[arctic$Year == i,], as.Date(paste(Year, Month, Day, sep="-")))
  days <- julian(date.list, origin = as.Date(ISOdate(i, 1, 1)) - 1)

  obj.lst <- c(obj.lst, list(list(extent = extent, count = count, days = days)))
}

names(obj.lst) <- paste("year", 1978:2014, sep="")

plot(x=c(1,366), y=c(1976, 2014), type="n")
for (i in 1:length(obj.lst)) {
  points(x = obj.lst[[i]]$days, y=rep(i + 1975, length(obj.lst[[i]]$days)), type = "o", pch=20)
}

# Remove years 1976, 85, 86
obj.lst[[1]] <- NULL
obj.lst[[9]] <- NULL
obj.lst[[9]] <- NULL

####### Create a loop to apply the above code to all years
require(fda)
NBASIS <- 300
NORDER <- 4

# One set of basis functions for all years
splinebasis <- create.bspline.basis(rangeval=c(1, 366),nbasis=NBASIS,norder=NORDER)
fdParobj<-fdPar(fdobj=splinebasis,Lfdobj=2,lambda=1)

predicted.mat <- matrix(NA, 34, 365)

for (i in 1:length(obj.lst)) {
  icefd <- smooth.basis(argvals=c(obj.lst[[i]]$days), y=obj.lst[[i]]$extent, fdParobj=fdParobj)$fd
  predicted.mat[i,] <- predict(icefd, newdata = 1:365)
}

par(mfrow = c(1,2))
plot(obj.lst$year1999$extent ~ obj.lst$year1999$days, type="n",
     ylim = c(min(arctic$Extent), max(arctic$Extent)))

for (i in 1:length(obj.lst)) {
  lines(obj.lst[[i]]$extent ~ obj.lst[[i]]$days, type="l", col = i)
}

plot(obj.lst$year1999$extent ~ obj.lst$year1999$days, type="n",
     ylim = c(min(arctic$Extent), max(arctic$Extent)))
for (i in 1:length(obj.lst))
  lines(predicted.mat[i,] ~ c(1:365), type="l", col = i)

Jan<-c(1,31)
Feb<-c(31,59)
Mar<-c(59,90)
Apr<-c(90,120)
May<-c(120,151)
Jun<-c(151,181)
Jul<-c(181,212)
Aug<-c(212,243)
Sep<-c(243,273)
Oct<-c(273,304)
Nov<-c(304,334)
Dec<-c(334,365)


intervals=rbind(Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)
abline(v=unique(as.vector(intervals)),lwd=2,col="grey")
text(x=mean(intervals[1,]),y=5,"Jan",cex=1.2)
text(x=mean(intervals[2,]),y=5,"Feb",cex=1.2)
text(x=mean(intervals[3,]),y=5,"Mar",cex=1.2)
text(x=mean(intervals[4,]),y=5,"Apr",cex=1.2)
text(x=mean(intervals[5,]),y=5,"May",cex=1.2)
text(x=mean(intervals[6,]),y=15,"Jun",cex=1.2)
text(x=mean(intervals[7,]),y=15,"Jul",cex=1.2)
text(x=mean(intervals[8,]),y=15,"Aug",cex=1.2)
text(x=mean(intervals[9,]),y=15,"Sep",cex=1.2)
text(x=mean(intervals[10,]),y=15,"Oct",cex=1.2)
text(x=mean(intervals[11,]),y=15,"Nov",cex=1.2)
text(x=mean(intervals[12,]),y=15,"Dec",cex=1.2)

## BEGIN: Rainbow object ##
x <- 1:365
y <- t(predicted.mat)
colnames(y) <- as.character(c(1977:1984, 1987:2014))
rownames(y) <- as.character(1:365)
time <- ts(start=1979, end = 2012)
xname <- "Year"
yname <- "Ice extent"
IceExtentRB <- list(x = x, y = y, time = time, xname = xname, yname = yname)
attr(IceExtentRB, "class") <- c("sfts", "fts", "fds")

str(IceExtentRB)
## END: Rainbow object ##

#Build common argval fd object from predicted.mat
y<-t(predicted.mat)
splinebasis <- create.bspline.basis(rangeval=c(1, 366),nbasis=NBASIS,norder=NORDER)
fdParobj<-fdPar(fdobj=splinebasis,Lfdobj=2,lambda=.000001)
yfd<-smooth.basis(argvals=1:365, y=y, fdParobj=fdParobj)

plot(yfd,ask=T)

require(fda.usc)

fd1<-fdata(yfd$fd)
km4<-kmeans.fd(fd1,ncl=4)

require(beanplot)
beanplot(as.numeric(yfd$fd$fdnames$reps)~km4$cluster,col="beige")

#require(devtools)
#install_github("vinhtantran/monoClust")
#require(monoClust)


# sourceDir <- function(path, trace = TRUE, ...) {
#   for (nm in list.files(path, pattern = "\\.R$")) {
#     if(trace) cat(nm,":")
#     source(paste(path, "/", nm, sep = ""), ...)
#     if(trace) cat("\n")
#   }
# }
#
# sourceDir("MonoClust/R", trace=F)
# source("MonoClust/R/print.MonoClust.R")
# source("MonoClust/R/labels.MonoClust.R")
# source("MonoClust/R/PCAmixcode.R")
# source("MonoClust/R/plots.R")
# source("MonoClust/R/text.MonoClust.R")
# source("MonoClust/R/plot.MonoClust.R")
# source("MonoClust/R/predict.MonoClust.R")
# source("MonoClust/R/MonoClust.R")
# source("MonoClust/R/cv.test.mse.R")
# source("MonoClust/R/cv.test.R")
# source("MonoClust/R/cv.plot.R")
# source("MonoClust/R/perm.test.R")

data1<-t(y)

MAXNODES <- 10

cp.table<-cv.test.mse(data.frame(data1),maxnodes = MAXNODES)



cv.plot(cp.table, main="")
mincp_order<-cp.table[,1]==min(cp.table[,1])
abline(h=(cp.table[mincp_order,1] + cp.table[mincp_order,2]), lty=2)

# Check if MonoClust can work on a list of selected variables
out <- MonoClust(data.frame(data1), variables = 16:35, nclusters=5)
tree.with.p <- perm.test(out,data.frame(data1))
plot(tree.with.p)

lwdid<-rep(1,34)
lwdid[out$centroids[,1]]<-3

plot(yfd,col=as.numeric(as.factor(out$Membership)),lty=as.numeric(as.factor(out$Membership)))
plot(yfd$fd[out$medoids],lwd=3,lty=1,col=as.numeric(as.factor(as.numeric(names(out$medoids)))),add=T)


require(beanplot)
beanplot(as.numeric(yfd$fd$fdnames$reps)~out$Membership,col="beige")

plot(yfd,col=as.numeric(as.numeric(yfd$fd$fdnames$reps)>2009)+1,lty=1)




x <- 1:365
y <- t(predicted.mat)
colnames(y) <- as.character(c(1979:1986, 1989:2014))
rownames(y) <- as.character(1:365)
time <- ts(start=1979, end = 2012)
xname <- "Year"
yname <- "Ice extent"
IceExtentRB <- list(x = x, y = y, time = time, xname = xname, yname = yname)
attr(IceExtentRB, "class") <- c("sfts", "fts", "fds")

str(IceExtentRB)


#########
data(ruspini)
ruspini <- cbind(ruspini, z = ruspini$y)
out <- MonoClust(ruspini, nclusters=4)
out <- MonoClust(ruspini, variables = 2:3, nclusters=4)
out <- MonoClust(ruspini, variables = c(1, 3), nclusters=4)
ruspini

out <- MonoClust(ruspini, variables = 2, nclusters=4)

#### 4/5
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.R$")) {
    if(trace) cat(nm,":")
    source(paste(path, "/", nm, sep = ""), ...)
    if(trace) cat("\n")
  }
}

sourceDir("../../R", trace=F)
source("../R/print.MonoClust.R")
source("../R/labels.MonoClust.R")
source("../R/PCAmixcode.R")
source("../R/plots.R")
source("../R/text.MonoClust.R")
source("../R/plot.MonoClust.R")
source("../R/predict.MonoClust.R")
source("../R/MonoClust.R")
source("../R/cv.test.mse.R")
source("../R/cv.test.R")
source("../R/cv.plot.R")
source("../R/perm.test.R")
