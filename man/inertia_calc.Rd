% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inertia_calc.R
\name{inertia_calc}
\alias{inertia_calc}
\title{Cluster Inertia Calculation}
\usage{
inertia_calc(x)
}
\arguments{
\item{x}{Distance matrix, not an object of some distance measure.}
}
\value{
Inertia value of the matrix, formula in Chavent (1998). If \code{x} is a
single number, return 0.
}
\description{
Calculate inertia for a given subset of the distance matrix from the original
data set provided to \code{x}. Assumes that distance matrices are stored as
matrices and not distance objects.
}
\examples{
data(iris)

# Euclidean distance on first 20 rows of the 4 continuous variables
dist_mat <- as.matrix(dist(iris[1:20, 1:4]))
inertia_calc(dist_mat)
}
