## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6
)

## ----ruspinivignette----------------------------------------------------------
library(monoClust)
library(cluster)
data(ruspini)
ruspini4c <- MonoClust(ruspini, nclusters = 4)
ruspini4c

## ----ruspinip, fig.height=5, fig.cap="Binary partitioning tree with three splits, four clusters for *ruspini* data."----
plot(ruspini4c)

## ----cptable------------------------------------------------------------------
set.seed(12345)
cp.table <- cv.test(ruspini, fold = 5, minnodes = 1, maxnodes = 10)
cp.table

## ----cvv, message=FALSE, fig.cap="The choice of clusters for Ruspini data made by 10-fold CV where *minCV* selects 10 clusters and *1SE* selects 4. The error bars are the $\\overline{MSE} \\pm 1SE$ and the choice of 4 clusters, the simplest solution within 1 standard error of the minimum error estimate (the dashed lines coincide with the bar at 10 clusters) is highlighted with a $\\times$."----
library(dplyr)
library(ggplot2)
ggcv(cp.table) +
  geom_hline(aes(yintercept = min(lower1SD)), color = "red", linetype = 2) +
  geom_hline(aes(yintercept = min(upper1SD)), color = "red", linetype = 2) +
  geom_point(aes(x = ncluster[4], y = MSE[4]), color = "red", size = 2) +
  geom_point(aes(x = ncluster[4], y = MSE[4]), color = "red", size = 5, shape = 4)

## ----hyptestv, fig.height=5, fig.width = 10, fig.cap="Binary partitioning tree with five splits, six clusters, but one split should be pruned based on its p-value of 0.8."----
ruspini6c <- MonoClust(ruspini, nclusters = 6)
ruspini6c.pvalue <- perm.test(ruspini6c, data = ruspini, method = "sw", rep = 1000)
plot(ruspini6c.pvalue, branch = 1, uniform = TRUE)

## ----sensit2008plot, fig.cap="Splitting rule for the four-cluster solution. The color at the node can be set by `cols` argument. They match the ones in Figure \\@ref(fig:PCPellipsev)."----
data(wind_sensit_2008)
# For the sake of speed in the example
wind_reduced_2008 <- wind_sensit_2008[sample.int(nrow(wind_sensit_2008), 50), ]
sensit042008 <- MonoClust(wind_reduced_2008, nclusters = 4, cir.var = 3)

## ----PCPellipsev, fig.cap = "PCP with the circular variable (*WDIR*) depicted as an ellipse. The geographical direction is noted and the ellipse is rotated to facilitate understanding of clusters."----
ggpcp(data = wind_reduced_2008, 
      circ.var = "WDIR",
      rotate = pi / 4 + 0.6,
      order.appear = c("WDIR", "has.sensit", "WS"),
      clustering = sensit042008$membership, 
      medoids = sensit042008$medoids,
      alpha = 0.5,
      cluster.col = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"),
      show.medoids = TRUE)

