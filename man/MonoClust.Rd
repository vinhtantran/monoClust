\name{MonoClust}
\alias{MonoClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MonoClust is the main function of the package.
}
\description{
MonoClust creates a MonoClust object after partitioning the data set using MonoClust clustering.
}
\usage{
MonoClust(toclust, distmethod = NULL, labels = as.character(1:length(toclust[, 1])),
          digits = options("digits")$digits, nclusters = nrow(toclust),
          minbucket = round(minsplit/3), minsplit = 5,
          corders = 2, alpha = 0.05, perm.test = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{toclust}{
the data set
}
  \item{distmethod}{
distance method to use with the data set. The default value is the Euclidean distance for quantitative variables and the ... distance for categorical variables.
}
  \item{labels}{
define names of variables.
}
  \item{digits}{
significant numbers shown in the result.
}
  \item{nclusters}{
number of clusters created.
}
  \item{minbucket}{
the minimum number of observations in any terminal <leaf> node. If only one of minbucket or minsplit is specified, the code either sets minsplit to minbucket*3 or minbucket to minsplit/3, as appropriate.
}
  \item{minsplit}{
the minimum number of observations that must exist in a node in order for a split to be attempted.
}
  \item{corders}{
%%     ~~Describe \code{corders} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{perm.test}{
%%     ~~Describe \code{perm.test} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (toclust, distmethod = NULL, labels = as.character(1:length(toclust[,
    1])), digits = options("digits")$digits, nclusters = nrow(toclust),
    minbucket = round(minsplit/3), minsplit = 5, corders = 2,
    alpha = 0.05, perm.test = FALSE)
{
    if (minbucket >= minsplit) {
        cat("minbucket must be less than minsplit")
        return(0)
    }
    assign(".MonoClustwarn", 0, envir = .GlobalEnv)
    weights <- rep(1, nrow(toclust))
    quali_ordered <- NULL
    i <- sapply(toclust, is.character)
    toclust[i] <- lapply(toclust[i], function(x) as.factor(x))
    i <- sapply(toclust, is.factor)
    if (sum(i) != 0 & is.null(distmethod)) {
        distmethod <- "gower"
    }
    if (sum(is.na(toclust))) {
        imputed <- mice(toclust)
        cat("\nData contain missing values mice() used for imputation")
        cat("\nSee mice() help page for more details")
        cat("\nMissing cells per column:")
        print(imputed$nmis)
        toclust <- complete(imputed)
    }
    toclust0 <- toclust
    factors <- sapply(toclust, is.factor)
    qualtog <- ifelse(sum(factors) > 0, 1, 0)
    quanttog <- ifelse(sum(!factors) > 0, 1, 0)
    extra <- sum(factors)
    extracols <- ncol(toclust) + 1:(extra * (corders - 1))
    if (qualtog) {
        numbyvar <- sapply(toclust[, factors], function(xxx) length(levels(xxx)))
        catnames <- colnames(toclust[, factors])
        catrepvarlevel <- unlist(sapply(toclust[, factors], function(x) sort(levels(x))))
        names(catrepvarlevel) <- rep(catnames, numbyvar)
        PCA <- PCAmix(X.quanti = toclust[, !factors], X.quali = toclust[,
            factors], graph = FALSE)
        catrepvarlevelordered <- character(0)
        cuts_quali <- numeric(0)
        position <- 0
        for (j in 1:corders) {
            catrepvarlevel <- catrepvarlevel[order(PCA$categ.coord[,
                j])]
            catrepvarlevelordered <- cbind(catrepvarlevelordered,
                catrepvarlevel)
            for (p in 1:length(catnames)) {
                position <- position + 1
                quali_ordered[[position]] <- catrepvarlevel[names(catrepvarlevel) ==
                  catnames[p]]
                fcolumn <- toclust[, which(factors)[p]]
                var_order <- quali_ordered[[position]]
                cuts_quali <- cbind(cuts_quali, sapply(fcolumn,
                  function(xxx) match(xxx, var_order) + 1))
            }
        }
        qual_quant <- cuts_quali - 1
        toclust[, factors] <- qual_quant[, 1:extra]
        toclust[, extracols] <- qual_quant[, -c(1:extra)]
    }
    quantis <- !factors
    quantis <- c(quantis, rep(FALSE, length(extracols)))
    findclosest <- function(col) {
        sapply(col, function(x) ifelse(x == max(col), x + 1,
            min(col[which(col - x > 0)])))
    }
    if (quanttog) {
        cuts_quant <- apply(toclust[, quantis], 2, findclosest)
    }
    cuts <- toclust
    if (qualtog) {
        cuts[, c(which(factors), extracols)] <- cuts_quali
    }
    if (quanttog) {
        cuts[, which(quantis)] <- cuts_quant
    }
    if (!qualtog) {
        catnames <- character(0)
    }
    labs <- labels
    nvars <- length(toclust[1, ])
    nobs <- length(labs)
    distmats <- matrix()
    data <- as.data.frame(toclust)
    distcols <- c(1:ncol(toclust0), rep(which(factors), corders -
        1))
    if (length(catnames)) {
        othercolnames <- paste(colnames(toclust)[rep(which(factors),
            corders - 1)], rep(c(2:corders), each = sum(factors)),
            sep = "*~*")
        currcolnames <- paste(colnames(toclust)[which(factors)],
            1, sep = "*~*")
        colnames(toclust)[which(factors)] <- currcolnames
        colnames(toclust)[-c(1:ncol(toclust0))] <- othercolnames
        catnames <- c(currcolnames, othercolnames)
    }
    colnames(cuts) <- colnames(toclust)
    distmat0 <- daisy(toclust0)
    distmats <- as.matrix(distmat0)
    members <- 1:nobs
    assign(".Cloc", rep(1, nobs), envir = .GlobalEnv)
    assign(".Cluster_frame", data.frame(number = 1, var = "<leaf>",
        n = nobs, wt = sum(weights[members]), inertia = inertia(distmats[members,
            members]), bipartvar = "NA", bipartsplitrow = NA,
        bipartsplitcol = NA, inertiadel = 0, yval = 1, medoid = med(members,
            distmats), category = NA, cut = NA, loc = 0.1, stringsAsFactors = FALSE,
        split.order = 0), envir = .GlobalEnv)
    split.order <- 1
    while (sum(.Cluster_frame$var == "<leaf>") < nclusters) {
        check <- checkem(toclust, cuts, distmats, catnames, weights,
            minsplit, minbucket, split.order)
        split.order <- split.order + 1
        if (check == 0) {
            break
        }
    }
    rownames(.Cluster_frame) <- .Cluster_frame$number
    .Cluster_frame2 <- .Cluster_frame[, -1]
    textfxn <- function(yval, dev, wt, ylevel, digits, n, meds,
        names, use.n) {
        paste("\n  n=", n, "\n           M=", meds, sep = "")
    }
    var <- .Cluster_frame2$var
    cattog <- .Cluster_frame2$category
    splits <- which(var != "<leaf>")
    cat_splits <- which(var != "<leaf>" & cattog == 1)
    labs <- getlevels(splits, cats = cat_splits, varnames = var,
        frame = .Cluster_frame2, catnames = catnames, quali_ordered = quali_ordered,
        digits = digits)
    colnames(.Cluster_frame2)[4] <- "dev"
    .Cluster_frame2 <- .Cluster_frame2[, c(1, 12, 2, 3, 4, 5,
        6, 7, 8, 9, 10, 11, 13, 14)]
    dendfxns <- list(text = textfxn)
    Terms <- colnames(toclust0)
    centroids <- find.centroid(toclust0, qualtog, quanttog)
    medoids <- .Cluster_frame2[.Cluster_frame2$var == "<leaf>",
        "medoid"]
    names(medoids) <- rownames(.Cluster_frame2[.Cluster_frame2$var ==
        "<leaf>", ])
    rpartobj <- list(frame = .Cluster_frame2, labels = labs,
        labelsnum = labs, functions = dendfxns, qualordered = quali_ordered,
        Membership = .Cloc, Dist = distmats, Catnames = catnames,
        terms = Terms, centroids = centroids, medoids = medoids)
    class(rpartobj) <- c("MonoClust", "rpart")
    rm(list = c(".Cluster_frame", ".Cloc"), envir = globalenv())
    return(rpartobj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
