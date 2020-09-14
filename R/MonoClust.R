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
#'   variable (`cir.var != NULL`).
#' @param labels Displayed names of observations. If not specified and data
#' frame has row names, row names will be used. Otherwise, it is row index.
#' @param digits Significant decimal number printed in the output.
#' @param nclusters Number of clusters created. Default is 2.
#' @param minsplit The minimum number of observations that must exist in a node
#'   in order for a split to be attempted. Default is 5.
#' @param minbucket The minimum number of observations in any terminal leaf
#'   node. If not specified, it is set to minsplit/3.#
#' @param perm.test Whether or not to make a permutation test as stopping
#'   criterion while clustering. Default is FALSE.
#' @param alpha Value applied specifically to permutation test. Only valid when
#'   `perm.test = TRUE`.
#' @param ran This parameter should never be used.
#'
#' @return MonoClust object, an extension of rpart object.
#' @importFrom dplyr `%>%`
#' @export
#'
#' @references
#' 1. Chavent, M. (1998). A monothetic clustering method. Pattern Recognition
#' Letters, 19(11), 989–996. https://doi.org/10.1016/S0167-8655(98)00087-7
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
#' # data with categorical variable
#' data(iris)
#' iris_cat <- MonoClust(iris, nclusters = 3)
#' iris_cat
#'
#' # data with circular variable
#' data(wind_sensit_2007)
#'
#' circular_wind <- MonoClust(wind_sensit_2007, cir.var = 3, nclusters = 2)
#' circular_wind
MonoClust <- function(toclust,
                      cir.var = NULL,
                      variables = NULL,
                      distmethod = NULL,
                      labels = NULL,
                      digits = options("digits")$digits,
                      nclusters = 2,
                      minsplit = 5,
                      minbucket = round(minsplit / 3),
                      perm.test = FALSE,
                      alpha = 0.05,
                      ran = 0) { # Tan 4/16 Added to control recursive call

  if (!is.data.frame(toclust)) {
    stop("toclust must be a data frame.")
  } else {
    if (is.null(labels)) {
        labels <- ifelse(tibble::has_rownames(toclust),
                         rownames(toclust),
                         as.character(seq_len(nrow(toclust))))
    }
    toclust <- dplyr::as_tibble(toclust)
  }

  ## MOVE: Tan, 12/14, move to the top to save some calculations if bad
  ## parameters are transferred
  ## Ensure that important options make sense
  if (minbucket >= minsplit) {
    stop("minbucket must be less than minsplit.")
  }

  # ADD: Tan, 9/9/2020, add to check that no categorical variables existed
  if (!all(purrr::map_lgl(toclust, is.numeric))) {
    stop("Function does not support categorical variables yet.")
  }

  ## Tan, 5/27/16, argument checking (variables)
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
    # Tan 4/16/17, if there is a circular variable in the data set, use Gower's
    ## distance unless otherwise specified.
    distmethod <- ifelse(!is.null(cir.var), "gower", "euclidean")
  } else {
    distmethod <- match.arg(distmethod, c("euclidean", "manhattan", "gower"))
  }

  ## Tan, 4/10/17, argument checking (circular variables)
  if (!is.null(cir.var)) {
    if (!is.vector(cir.var)) {
      stop("Circular variables need to be a vector of variable names or
           indices.")
    }
    toclust1 <- toclust[cir.var]
    if (!is.data.frame(toclust1)) {
      stop("Undefined variables selected.")
    }
  }

  # R#EMOVE: Tan, 9/8/20. Don't see the benefits of this.
  ## Switch so we only give off one warning.
  # assign(".MonoClustwarn", 0, envir = .GlobalEnv)
  ## Right now, each observation has equal weight. This could be made into an
  ## option.
  weights <- rep(1, nrow(toclust))

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  ## Categorical Variable Ordering
  # quali_ordered <- NULL

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  ## REMOVE: Tan, 4/6/15. These codes haven't show usefulness, but R v1.3.1
  ## changed output of lapply to list, fail running.
  ## ADD: Tan, 9/6/20. Change to dplyr.
  ## See which variables are character and convert them to factors.
  # toclust <- purrr::modify_if(toclust, is.character, as.factor)

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  ## If we have categorical variables we want to use Gower's distance unless
  ## otherwise specified.
  # if (any(purrr::map_lgl(toclust, is.factor)) && is.null(distmethod)) {
  #   distmethod <- "gower"
  # }

  ## For now, impute missing values so we have no missing data
  ## Emit warning that values were imputed.
  if (any(is.na(toclust))) {
    imputed <- mice::mice(toclust)
    cat("\nData contain missing values mice() used for imputation")
    cat("\nSee mice() help page for more details")
    cat("\nMissing cells per column:")
    print(imputed$nmis)
    toclust <- mice::complete(imputed)
  }

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  ## Determine whether we need to use PCAmix, or whether we only have
  ## quantitative variables.
  # factors <- purrr::map_lgl(toclust, is.factor)
  # qualtog <- any(factors)
  # quanttog <- any(!factors)

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  # If has at least one categorical variable
  # extracols <- NULL
  # catnames <- character(0)
  # TODO: Candidate to move to a function
  # Input: toclust, factors, corders
  # Output: toclust, quantis
  # if (qualtog) {
  #   ## Since we are going to have multiple factor orderings,
  #   ## we need to consider each of these orderings a separate "explanatory"
  #   ## variable, and determine the number of extra columns we will need based on
  #   ## the number of orderings and the number of factors in the dataset.
  #   num_factors <- sum(factors)
  #   extracols <- ncol(toclust) + 1:(num_factors * (corders - 1))
  #
  #   ## Use PCAmix with quantitative variables the non-factors and qualitative
  #   ## variables the factors.
  #   pca <- PCAmixdata::PCAmix(X.quanti = toclust[, !factors],
  #                             X.quali = toclust[, factors],
  #                             graph = FALSE)
  #
  #   ## Get the number of levels of each factor.
  #   numbyvar <- purrr::map_int(toclust[, factors], nlevels)
  #
  #   ## Get each factor variable name.
  #   catnames <- colnames(toclust)[factors]
  #
  #   ## Sort the levels alphabetically as this is the way they will come out of
  #   ## PCAmix.
  #   catrepvarlevel <- purrr::flatten_chr(purrr::map(toclust[, factors],
  #                                                   ~ sort(levels(.x))))
  #   names(catrepvarlevel) <- rep(catnames, numbyvar)
  #
  #   ## Set up variables to be filled in the loop
  #   ## This method seems a bit odd, but PCAmix looks at all levels of all
  #   ## variables simultaneously, and some work must be done to keep the levels
  #   ## with the factor they belong to.
  #   catrepvarlevelordered <- character(0)
  #   cuts_quali <- numeric(0)
  #   position <- 0
  #
  #   for (j in 1:corders) {
  #     ## Do this for each dimension that we want 1-5 (number of orderings)
  #     ## Order and put them in columns
  #     catrepvarlevel <- catrepvarlevel[order(pca$categ.coord[, j])]
  #     catrepvarlevelordered <- cbind(catrepvarlevelordered, catrepvarlevel)
  #
  #     for (p in seq_len(num_factors)) {
  #       position <- position + 1
  #       quali_ordered[[position]] <-
  #         catrepvarlevel[names(catrepvarlevel) == catnames[p]]
  #
  #       ## We want the orders and the cuts
  #       ## The cuts are the ordered levels that we can use to partition the
  #       ## dataset. They are numeric, but correspond to a division between two
  #       ## levels of a categorical variable according to a particular ordering
  #       ## from PCAmix.
  #       fcolumn <- toclust[, which(factors)[p]]
  #       cuts_quali <-
  #         cbind(cuts_quali,
  #               purrr::map_dbl(fcolumn,
  #                              ~ match(.x, quali_ordered[[position]]) + 1))
  #     }
  #   }
  #
  #   qual_quant <- cuts_quali - 1
  #   toclust[, factors] <- qual_quant[, 1:num_factors]
  #   toclust[, extracols] <- qual_quant[, -c(1:num_factors)]
  #
  # }
  ## Now, we have extra columns, we have to make sure that we do not consider
  ## them quantitative.
  # quantis <- !factors
  # quantis <- c(quantis, rep(FALSE, length(extracols)))

  # CLUSTERING ON CIRCULAR VARIABLE
  # Keep an original copy of toclust before applying circular sp;it
  toclust0 <- toclust

  if (!is.null(cir.var)) {

    # This only works on **one** circular variable
    next_value <- find_closest(toclust[, cir.var])
    bestcircsplit <- list(hour = 0, minute = 0, inertia = 0)
    ## Tan 4/16/17, discreetly find the first split for the circular variable to
    ## nail the the first split
    if (ran == 0) {
      variable <- toclust[, cir.var]
      min_value <- min(variable)
      min_inertia <- 9999
      # Check all starting value to find the best split
      # The best split will be recorded in output with the pivot is hour
      while (min_value < max(variable)) {
        variable_shift <- cshift(variable, -min_value)
        # TODO: Have to call MonoClust with ran = 1 to find the best split for
        # the shifted variables. Should be improved.
        out <- MonoClust(data.frame(variable_shift),
                         cir.var = 1,
                         nclusters = 2,
                         ran = 1)
        cut <- out$frame$cut[1]
        inertia <-
          inertia_calc(out$dist[which(out$Membership == 2),
                                which(out$Membership == 2)]) +
          inertia_calc(out$dist[which(out$Membership == 3),
                                which(out$Membership == 3)])

        if (min_inertia > inertia) {
          min_inertia <- inertia
          bestcircsplit <- list(hour = min_value,
                                minute = ifelse((cut + min_value) < 360,
                                                cut + min_value,
                                                cut + min_value - 360),
                                intertia = inertia)
        }

        # Increase min_value to the next higher value
        min_value <- as.numeric(next_value[which(variable == min_value)])[1]
      }

      # Shift the circular variable to the hour pivot. That turns the circular
      # to linear variable. The first best split would be the minute split,
      # which was already found, together with hour, to be the best arcs. Will
      # shift back later by modifying .Cluster_frame2
      toclust[, cir.var] <- cshift(variable, -bestcircsplit$hour)
    }
  }

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  # Clustering on common quantitative variables
  # cuts <- toclust
  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  # if (qualtog) cuts[, c(which(factors), extracols)] <- cuts_quali

  # MODIFY: Tan, 9/9/20. Remove categorical variable for now.
  cuts <- purrr::map_dfc(toclust, find_closest)
  # cuts[, which(quantis)] <- cuts_quant

  # if (!qualtog) catnames <- character(0)

  ## Variables that are simple derivatives of inputs that will be used a lot.
  # labs <- labels
  nobs <- nrow(toclust)
  # nvars <- dim(toclust)[2]

  # dismat <- matrix()
  # data <- as.data.frame(toclust)

  ## which column to use in the distance matrix
  # distcols <- c(1:ncol(toclust0),rep(which(factors),corders-1))

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  ## Since we have multiple categorical orderings for each factor variable,
  ## We need to label each of these orderings. I chose to do it in a kind of
  ## bizarre way.
  ## The name followed by *~* and then the number.
  ## This is odd so that we can replace the *~* in the end which is a unique
  ## character combination and likely wont interfere with any of the variable
  ## names.
  # if (length(catnames) > 0) {
  #   othercolnames <- paste(colnames(toclust)[rep(which(factors), corders - 1)],
  #                          rep(c(2:corders), each = sum(factors)),
  #                          sep="*~*")
  #   currcolnames  <- paste(colnames(toclust)[which(factors)], 1, sep="*~*")
  #   colnames(toclust)[which(factors)] <- currcolnames
  #   colnames(toclust)[-c(1:ncol(toclust0))] <- othercolnames
  #   catnames <- c(currcolnames, othercolnames)
  # }

  # REMOVE: Tan 9/8/20, no need to assign colnames due to tibble's benefits
  # colnames(cuts) <- colnames(toclust)

  ## Tan 4/10/17, add metric argument to daisy and circular distance
  if (!is.null(cir.var)) {
    # I have to split because R coerce the distance matrix to one value when
    # using ifelse with a 0.
    # distmat1: distance matrix of pure quantitative variables
    if (ncol(toclust0[, -cir.var]) != 0) {
      distmat1 <- cluster::daisy(toclust0[, -cir.var], metric = distmethod) *
        dim(toclust0[, -cir.var])[2]
    } else {
      distmat1 <- 0
    }
    # distmat2: distance matrix of circular variables
    # TODO cir.var has to be one value
    distmat2 <- circ_dist(toclust0[, cir.var])

    # distmat0: combined from 1 and 2 as the mean distances
    distmat0 <- (distmat1 + distmat2) / dim(toclust0)[2]
  } else
    distmat0 <- cluster::daisy(toclust0, metric = distmethod)

  dismat <- as.matrix(distmat0)

  members <- seq_len(nobs)

  ## Set up a vector containing each observation's membership. Put into global
  ## environment, but this will be deleted at the end of this function. Using
  ## global environment allows us to modify things recursively as we partition
  ## clusters.
  assign(".Cloc", rep(1, nobs), envir = .GlobalEnv)

  ## Likewise, set up the first (entire dataset) cluster in our Cluster frame
  ## where we keep track of each of the clusters and the partitioning.
  assign(".Cluster_frame",
         data.frame(number = 1,
                    var = "<leaf>",
                    n = nobs,
                    wt = sum(weights[members]),
                    inertia = inertia_calc(dismat[members, members]),
                    bipartvar = "NA",
                    bipartsplit_row = NA,
                    bipartsplitcol = NA,
                    inertiadel = 0,
                    yval = 1,
                    medoid = medoid(members, dismat),
                    category = NA,
                    cut = NA,
                    loc = 0.1,
                    stringsAsFactors = FALSE,
                    split.order = 0),
         envir = .GlobalEnv)

  split_order <- 1
  done_running <- FALSE
  ## This loop runs until we have nclusters, have exhausted our observations or
  ## run into our minbucket/minsplit restrictions.
  while (sum(.Cluster_frame$var == "<leaf>") < nclusters & !done_running) {

    # MODIFY: Tan, 9/9/20. Remove categorical variable for now.
    checkem_ret <- checkem(toclust, cuts, .Cluster_frame, .Cloc, dismat,
                           variables, weights, minsplit, minbucket, split_order)
    split_order <- split_order + 1

    # Use cloc because it is only ran in splitter
    if (!identical(.Cloc, checkem_ret$cloc)) {
      .Cluster_frame <- checkem_ret$frame
    } else {
      done_running <- TRUE
    }

    .Cloc <- checkem_ret$cloc
  }

  ## Most of the rest of the function does some bizarre text operations
  ## the reason for this is because I stole a lot of code from rpart,
  ## so we need to follow their text and labeling conventions to that our
  ## objects which inherit from rpart can print and plot correctly.

  tibble::column_to_rownames()

  ## Change the number column to rownames...
  rownames(.Cluster_frame) <- .Cluster_frame[[number]]
  .Cluster_frame[[number]] <- NULL

  ## This is what will print at each terminal node on the dendrogram
  ## (See plot.MonoClust).
  textfxn<-function(yval,dev,wt,ylevel,digits,n,meds,names, use.n){
    paste("\n  n=", n,"\n  M=",meds, sep="")
  }

  ## Seperate categorical and quantitative splits as the text and plot
  ## functions must treat them a bit differently
  var <-.Cluster_frame2$var
  cattog <- .Cluster_frame2$category

  splits<-which(var != "<leaf>")
  cat_splits<-which(var != "<leaf>" & cattog == 1)

  ## Piece together a vector of labels to be printed. Kind of a weird way to do this, but
  ## again, following rparts conventions, and we want to allow the user to have options
  ## regarding how to print inequalities.

  ## REMOVE: Tan, 12/14. The following lines are obviously useless. They have no use anywhere.
  # ineq<-rep(c('<','>='),length(splits))
  # level<-.Cluster_frame2$cut[splits]
  # level<-rep(level,each=2)
  # vars<-rep(var[splits],each=2)
  # labsnum <- c('root',paste(vars,ineq,level,sep=' '))

  ## Tan, 4/16/17
  ## Modify the Cluster_frame to shift the circular variable's cut back to the original values
  .Cluster_frame2[which(.Cluster_frame2$var == colnames(toclust)[cir.var]), "cut"] <-
    cshift(.Cluster_frame2[which(.Cluster_frame2$var == colnames(toclust)[cir.var]), "cut"],  bestcircsplit$hour)

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  ## MODIFY: Tan, 12/14. Change input of getlevels. If getlevels doesn't see the whole structure of output, we can't
  ## set correct left and right node labels. See new getlevels function for more details.
  # labs<-c('root',sapply(splits,getlevels,cats = cat_splits,varnames=var, frame=.Cluster_frame2,catnames=catnames,quali_ordered=quali_ordered))
  labs<-getlevels(splits, cats = cat_splits, varnames=var, frame=.Cluster_frame2, digits=digits)

  ## name a column what I probably should hav already named it, but I don't want to change all the code.
  colnames(.Cluster_frame2)[4] <-"dev"

  ## Reorder the columns so they print out nicely, again because I don't want to go back and change things.
  .Cluster_frame2 <- .Cluster_frame2[,c(1,12,2,3,4,5,6,7,8,9,10,11,13,14)]

  ## This follows somewhat odd rpart conventions.
  dendfxns<-list("text"=textfxn)

  ## Take variables' names out of original data set
  Terms <- colnames(toclust0)

  # MODIFY: Tan, 9/9/20. Remove categorical variable for now.
  ## ADD: Tan, 12/15, calculate the mean of each cluster
  centroids <- find.centroid(toclust0)

  ## ADD: Tan, 4/22/16, add medoids of each cluster
  medoids <- .Cluster_frame2[.Cluster_frame2$var =="<leaf>","medoid"]
  names(medoids) <- rownames(.Cluster_frame2[.Cluster_frame2$var =="<leaf>",])

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  # We will return a MonoClust object that also inherits from rpart with all of the neccesary components.
  ## MODIFY: Tan, 12/14, I don't know the difference between labels and labelsnum output
  ## although both of them are used in labels.MonoClust. Maybe for categorical variables?
  rpartobj<-list("frame"=.Cluster_frame2,"labels"=labs,"labelsnum" = labs, "functions"=dendfxns, Membership =.Cloc, Dist=dismat,
                 terms = Terms, # 12/9/14. Tan: add terms to keep track of variables name, in order to check the new data set
                 centroids = centroids, # 12/15. Tan: add centroids info, for prediction of quantitative
                 medoids = medoids, # 4/22/15. Tan: add medoids info
                 circularroot = list(var = cir.var, cut = bestcircsplit$hour) # 4/16/17. Tan: add the starting cut for circular variable
  )
  class(rpartobj)<-c("MonoClust","rpart")

  ## Get rid of our global assignments.
  rm(list = c(".Cluster_frame", ".Cloc"), envir = globalenv())

  return(rpartobj)

}

# MODIFY: Tan, 9/9/20. Remove categorical variable for now.
## ADD, Tan, 12/15, function to calculate the mean of each cluster.
## Currently do not work for categorical variables
find.centroid <- function(toclust) {

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  # Don't calculate if there is qualitative variable
  # if (qualtog) NA

  leaf <- .Cluster_frame[.Cluster_frame$var == "<leaf>", "number"]
  centroid.list <- as.numeric()
  centroid.list <- vector("list", length(leaf))
  for (i in seq_along(leaf)) {
    cluster <- as.matrix(toclust[.Cloc == leaf[i],])
    centroid <- apply(cluster, 2, mean)
    centroid.list[[i]] <- c(leaf[i], centroid)
  }
  centroid.list <- data.frame(matrix(unlist(centroid.list), nrow=length(leaf), byrow = T))
  colnames(centroid.list) <- c("cname", rep("", ncol(centroid.list)-1))
  centroid.list
}


abbreviate.t <- function(string,abbrev){
  ## This function abbreviates the text. We want to handle categories, which rpart is not really prepared for
  ## so we use somewhat uneccesary regular expressions to get rid of unwanted characters in our categorical orderings.
  charvec <- strsplit(string," ",fixed=TRUE)[[1]]
  abbreviated <- sapply(charvec,function(x)if( x != "<" & x != ">=" & grepl("^\\s*[^0-9]",x,perl=TRUE) ) substr(x,1,abbrev)
                        else x)
  paste(abbreviated,collapse=" ")
}

## MODIFY: Tan, 12/14. Change input as a vector instead of one value.
## getlevels needs to see the whole structure to decide correct labels
## NOTE: quantitative case has been changed significantly, categorical case is only slightly changed to refect the input change

# getlevels <- function(sind,cats,varnames,frame,catnames,quali_ordered){
#     ## A bit of a pain in the ass to get categorical ordering levels to print correctly.
#     ## To be honest, I forgot what the last part here does, but I am certain it is neccesary.
#     name <- varnames[sind]
#     if(sind %in% cats == 0){
#         level <- frame$cut[sind]
#         labs <- c(paste(name,"<",level,sep=" "),paste(name,">=",level,sep=" "))
#         return(labs)
#     } else{
#         qualind <- which(catnames==varnames[sind])[1]
#         toret <- c(paste(quali_ordered[[qualind]][1:(frame$cut[sind]-1)],collapse=" "),paste(quali_ordered[[qualind]][-c(1:(frame$cut[sind]-1))],collapse=" "))
#         return(toret)
#     }
# }

getlevels <- function(ind,cats,varnames,frame, digits=getOption('digits')){
  ## A bit of a pain in the ass to get categorical ordering levels to print correctly.
  ## To be honest, I forgot what the last part here does, but I am certain it is neccesary.

  # These codes are modified version of rpart:::labels.rpart
  lsplit <- rsplit <- character(length(ind))
  # If there exists quantitative cutpoint
  if (any(!ind %in% cats)) {
    sind <- ind[ind %in% cats == 0]
    name <- varnames[sind]
    level <- frame$cut[sind]
    lsplit[ind %in% cats == 0] <- paste(name,"<",round(level, digits),sep=" ")
    rsplit[ind %in% cats == 0] <- paste(name,">=",round(level, digits),sep=" ")
  }
  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  # If there exists categorical cutpoint
  # if (any(ind %in% cats)) {
  #   sind <- ind[ind %in% cats == 1]
  #   for (i in sind) {
  #     name <- varnames[i]
  #     qualind <- which(catnames==varnames[i])
  #     lsplit[which(ind == i)] <- paste(quali_ordered[[qualind]][1:(frame$cut[i]-1)],collapse=" ")
  #     rsplit[which(ind == i)] <- paste(quali_ordered[[qualind]][-c(1:(frame$cut[i]-1))],collapse=" ")
  #   }
  #
  # }

  node <- as.numeric(row.names(frame))
  parent <- match(node%/%2, node[ind])
  odd <- (as.logical(node%%2))

  labels <- character(nrow(frame))
  labels[odd] <- rsplit[parent[odd]]
  labels[!odd] <- lsplit[parent[!odd]]
  labels[1] <- "root"
  labels
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
splitter <- function(data, cuts, split_row, frame, cloc, dist, weights,
                     split_order = 0) {
  ## This function does the actual act of partitioning, given the row that is
  ## to be split "split_row"

  node_number <- frame$number[split_row]
  mems <- which(cloc == node_number)
  split <- c(frame$bipartsplit_row[split_row],
              frame$bipartsplitcol[split_row])

  # Extract data and cuts to split
  datamems <- data[mems, ]
  cutsmems <- cuts[mems, ]

  # Split into data rows of lower half (A) and upper half (B)
  mems_A <- mems[which(datamems[[split[2]]] < cutsmems[[split[1], split[2]]])]
  mems_B <- setdiff(mems, mems_A)

  # Tan, 2/17, call the permutation test function to test on the newly created
  # group perm.test condition should be uncommented later when successfully
  # tested
  # if (perm.test) {
  # ptest.result <- permtest(split[2], data, mems_A, mems_B)
  # ptest.result$aov.tab[1,6]
  # }

  # UPDATE: Tan, 9/13/20, cutpoint is simply value between two jump values
  # Tan, 10/3, cutpoint is the middle point between two closest points in two
  # clusters
  # datamems_A <- datamems[datamems[, split[2]] <
  #                         cutsmems[split[1], split[2]], split[2]]
  # datamems_B <- setdiff(datamems[, split[2]], datamems_A)
  # mid_cutpoint <- (max(datamems_A) + min(datamems_B)) / 2
  mid_cutpoint <- mean(datamems[[split[1], split[2]]],
                       cutsmems[[split[1], split[2]]])

  ## Make the new clusters.
  node_number_A <- node_number * 2
  node_number_B <- node_number * 2 + 1

  cloc[mems_A] <- node_number_A
  cloc[mems_B] <- node_number_B

  variable_name <- colnames(data)[frame$bipartsplitcol[split_row]]

  # REMOVE: Tan, 9/13/20. Remove categorical variable for now.
  ## This seperates the categorical variable from the level.
  ## Probably bad coding, parsing strings over and over.
  # if (grepl(variable_name, "*~*", fixed = TRUE)) {
  #   variable_name <- strsplit(variable_name, "*~*", fixed = TRUE)[[1]]
  # }

  # REMOVE: Tan, 9/9/20. Remove categorical variable for now.
  ## Is the split categorical?
  # if(variable_name %in% catnames) {
  #   frame[split_row,12] <<- 1
  # }else { frame[split_row,12] <<- 0 }
  frame$category[split_row]    <- 0

  ## The old cluster now changes some attributes after splitting.
  frame$var[split_row]         <- variable_name
  frame$bipartvar[split_row]   <- variable_name
  frame$cut[split_row]         <- mid_cutpoint # Use new cutpoint
  frame$split.order[split_row] <- split_order

  ## Tan, 12/14, trying to change the order of rows in Cluster_frame, so that
  ## the children will stay right after the parent, instead of at the very end
  ## of the table
  ## Reason: mimic rpart, help print, and plot.
  ## All the comment or new one will be noted carefully so that it can be
  ## rollbacked easily

  # Initialize two new empty rows
  new_frame <- frame[rep(1, 2), ]
  new_frame[, ] <- NA

  ## New cluster 1 gets some new attributes
  new_frame$number[1]  <- node_number_A
  new_frame$var[1]     <- "<leaf>"
  new_frame$n[1]       <- length(mems_A)
  new_frame$wt[1]      <- sum(weights[mems_A])
  new_frame$inertia[1] <- inertia_calc(dist[mems_A, mems_A])
  new_frame$yval[1]    <- 1 - frame[split_row, 9] / frame[1, 5]
  new_frame$medoid[1]  <- medoid(mems_A, dist)
  new_frame$loc[1]     <- frame[split_row, 14] - 1/nrow(frame)

  ## As does new cluster 2.
  new_frame$number[2]  <- node_number_B
  new_frame$var[2]     <- "<leaf>"
  new_frame$n[2]       <- length(mems_B)
  new_frame$wt[2]      <- sum(weights[mems_B])
  new_frame$inertia[2] <- inertia_calc(dist[mems_B, mems_B])
  new_frame$yval[2]    <- 1 - frame[split_row, 9] / frame[1, 5]
  new_frame$medoid[2]  <- medoid(mems_B, dist)
  new_frame$loc[2]     <- frame[split_row, 14] + 1/nrow(frame)

  # Insert two new rows right after split row
  frame <- dplyr::add_row(frame, new_frame, .after = split_row)

  return(list(frame = frame, cloc = cloc))
}

#' Find the Best Split
#'
#' Find the best split in terms of reduction in inertia for the transferred
#' node, indicate by row. We are keeping all of the information regarding which
#' clusters we have in. At this point, we want to find the terminal node with
#' the greatest change in inertia and bi-partition it. The way this is done is
#' incredibly inefficient. We check every single possible split. This should be
#' done much better since there is a single local maximum with respect to
#' inertia so we should use this property to not do an exhaustive check, but to
#' check the center the 25% and 75% and then search more selectively from there.
#' This is a fairly simple discrete optimization problem and could reduce
#' computation time .
#'
#' @param frame_row One row of the split tree as data frame.
#' @inheritParams checkem
#'
#' @return This function changes the frame in global environment, it's
#'   not supposed to return anything.
#' @importFrom foreach `%dopar%`
#' @keywords internal
find_split <- function(data, cuts, frame_row, cloc, dist, variables, minsplit,
                       minbucket) {

  node_number <- frame_row$number
  mems <- which(cloc == node_number)
  inertiap <- frame_row$inertia

  if (inertiap == 0 | frame_row$n < minsplit | frame_row$n == 1) {
    # Tan 9/24 This is one obs cluster. Set bipartsplit_row value 0 to stop
    # checkem forever.
    # MG, 9/25 I think this means we won't explore this node again. But make it
    # doesn't stop the search into other nodes.
    # Yup, but we waste resources by keep checking them again and again. The
    # "candidates" of checkem keeps getting longer (at most n) instead of just 2
    # new splits.
    frame_row$bipartsplit_row <- 0
    return(frame_row)
  }

  ## Subset the data and cut matricies
  # MODIFY: Tan, 7/3/16, add search space limit
  # datamems<-data[mems,]
  # cutsmems<-cuts[mems,]
  datamems <- data[mems, variables]
  cutsmems <- cuts[mems, variables]

  ## For each possible cut, calculate the inertia. This is where using a
  ## discrete optimization algorithm would help a lot.
  bycol <- foreach::foreach(index = seq_len(ncol(datamems)),
                   .combine = cbind,
                   .export = c("datamems", "cutsmems", "dist", "mems")) %dopar%
    {
      data_col <- dplyr::pull(datamems, index)
      cuts_col <- dplyr::pull(cutsmems, index)

      new_inertia <- purrr::map_dbl(cuts_col, function(x) {
        mems_A <- mems[which(data_col < x)]
        mems_B <- setdiff(mems, mems_A)
        ifelse(length(mems_A) * length(mems_B) == 0,
               NA,
               inertia_calc(dist[mems_A, mems_A]) +
                 inertia_calc(dist[mems_B, mems_B]))
      })

      return(new_inertia)
    }

  # for (i in seq_len(ncol(datamems))) {
  #
  #   data_col <- datamems[, i]
  #   cuts_col <- cutsmems[, i]
  #
  #   bycol <- cbind(bycol, sapply(cuts_col, function(x) {
  #     mems_A <- mems[which(data_col < x)]
  #     mems_B <- setdiff(mems, mems_A)
  #     ifelse(length(mems_A) * length(mems_B) == 0,
  #            NA,
  #            inertia_calc(dist[mems_A, mems_A]) + inertia_calc(dist[mems_B, mems_B]))
  #     }))
  # }
  # Difference between current cluster and the possible splits
  vals <- inertiap - bycol
  ## Say no difference if we have NA or infinite (happens when no split is possible)
  vals[!is.finite(vals) | is.na(vals)] <- 0

  ## This is the best split.
  maxval <- max(vals)

  ## This is the maximum inertia change index
  ind <- which((inertiap - bycol) == maxval, arr.ind = TRUE)

  ## Add one more column to check minbucket
  ind_1 <- cbind(ind, minbucket = TRUE)

  for (i in 1:nrow(ind_1)) {
    split <- ind_1[i, ]
    left_size <- sum(datamems[, split[2]] < cutsmems[[split[1], split[2]]])
    right_size <- length(mems) - left_size
    if (left_size < minbucket | right_size < minbucket) ind_1[i, 3] <- FALSE
  }

  ## Remove all row doesn't satisfy minbucket and make sure output is always a matrix even though it has only one row
  ind_1 <- matrix(ind_1[!ind_1[,3] == FALSE, ], ncol = 3)

  ## If multiple splits produce the same inertia change output a warning.
  #if(nrow(ind) > 1 & .MonoClustwarn==0){.MonoClustwarn <<- 1; warning("One or more of the splits chosen had an alternative split that reduced deviance by the same amount.")}

  # If there is at least one row that satisfies minbucket, pick the first one
  if (nrow(ind_1) != 0) {
    split <- ind_1[1, ]

    mems_A <- mems[which(datamems[, split[2]] < cutsmems[[split[1], split[2]]])]
    mems_B <- setdiff(mems, mems_A)

    # calculate our change in inertia
    inertiadel <- inertiap -
      inertia_calc(dist[mems_A, mems_A]) -
      inertia_calc(dist[mems_B, mems_B])

    ## Update our frame
    frame_row$bipartsplit_row <- split[1]
    # Save colname, not colindex
    frame_row$bipartsplitcol <- variables[split[2]]
    frame_row$inertiadel <- inertiadel
  } else  # Otherwise, stop as a leaf
    frame_row$bipartsplit_row <- 0

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
#' @param weights Weights on each observation. Hasn't been implemented in
#' exported function yet. Vector of 1 for all observations.
#' @param minsplit The minimum number of observations that must exist in a node
#'   in order for a split to be attempted.
#' @param split_order The control argument to see how many split has been done.
#' @inheritParams MonoClust
#'
#' @return It is not supposed to return anything because global environment was
#'   used. However, if there is nothing left to split, it returns 0 to tell the
#'   caller to stop running the loop.
#' @keywords internal
checkem <- function(data, cuts, frame, cloc, dist, variables, weights, minsplit,
                    minbucket, split_order) {

  ## Current terminal nodes
  candidates <- which(frame$var == "<leaf>" &&
                        is.na(frame$bipartsplit_row))
  ## Split the best one. Return to nada which never gets output.
  frame[candidates, ] <-
    purrr::map_dfr(candidates,
                   ~ find_split(data, cuts, frame[.x, ], cloc, dist, variables,
                                minsplit, minbucket))


  ## See which ones are left.
  candidates2 <- which(frame$var == "<leaf>" && frame$bipartsplit_row != 0)
  ## If nothing's left, stop running.
  check <- ifelse(length(candidates2) == 0, 0, 1)

  # otherwise, frame and cloc are not updated, cloc is used to check if done
  # running
  if (length(candidates2) != 0) {
    ## Find the best inertia change of all that are possible
    split_row <- candidates2[which.max(frame$inertiadel[candidates2])]

    ## Make new clusters from that cluster
    # splitter(split_row, data, cuts,dist,catnames,weights)
    splitter_ret <- splitter(data, cuts, split_row, frame, cloc, dist, weights,
                             split_order)
    # Tan, 9/24, in case there are more than one node equal to max
    # MG, 9/25, I thought the earlier code would make sure only one is identified
    # as top but that might not be true. It never caused a problem before.
    # Maybe because we fixed inertia fn, equal inertias occurred for small
    # clusters

    frame <- splitter_ret$frame
    cloc <- splitter_ret$cloc
  }


  return(list(frame = frame, cloc = cloc))
}
