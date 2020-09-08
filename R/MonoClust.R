#' Monothetic Clustering
#'
#' Creates a MonoClust object after partitioning the data set using Monothetic
#' Clustering.
#'
#' @param toclust Data set.
#' @param cir.var Index or name of the circular variable in the data set.
#' @param variables List of variables selected for clustering procedure. It
#'   could be a vector of variable indexes, or a vector of variable names.
#' @param distmethod Distance method to use with the data set. The default value
#'   is the Euclidean distance for quantitative variables and the Gower's
#'   distance for categorical variables.
#' @param labels Displayed names of variables.
#' @param digits Significant number of clusters created
#' @param nclusters Number of clusters created
#' @param minbucket The minimum number of observations in any terminal leaf
#'   node. If only one of minbucket or minsplit is specified, the code either
#'   sets minsplit to minbucket*3 or minbucket to minsplit/3, as appropriate.
#' @param minsplit The minimum number of observations that must exist in a node
#'   in order for a split to be attempted.
#' @param corders For data set with categorical variables, the number of
#'   dimensions should be used from PCAMix. Taken a value from 1 to 5.
#' @param alpha Value applied specifically to permutation test.
#' @param perm.test Whether or not to make a permutation test as stopping
#'   criterion while clustering.
#' @param ran This parameter should never be used
#'
#' @return MonoClust object, an extension of rpart object
#' @importFrom dplyr `%>%`
#' @export
#'
#' @examples
#' # Very simple data set
#' data(cluster::ruspini)
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
                      labels = as.character(seq_len(length(toclust[, 1]))),
                      digits = options("digits")$digits,
                      nclusters = nrow(toclust),
                      minbucket = round(minsplit / 3),
                      minsplit = 5,
                      corders = 2,
                      alpha = 0.05,
                      perm.test = FALSE,
                      ran = 0) { # Tan 4/16 Added to control recursive call

  ## MOVE: Tan, 12/14, move to the top to save some calculations if bad
  ## parameters are transferred
  ## Ensure that important options make sense
  if (minbucket >= minsplit) {
    stop("minbucket must be less than minsplit.")
  }

  ## Tan, 5/27/16, argument checking (variables)
  if (!is.null(variables)) {
    if (!is.vector(variables)) {
      stop("variables need to be a vector of variable names or indices.")
    }
    toclust1 <- toclust[variables]
    if (!is.data.frame(toclust1)) {
      stop("undefined variables selected.")
    }
    if (variables %in% colnames(toclust)) {
      variables <- which(colnames(toclust) == variables)
    }
  } else {
    variables <- seq_len(ncol(toclust))
  }

  ## Tan, 4/10/17, argument checking (circular variables)
  if (!is.null(cir.var)) {
    if (!is.vector(cir.var)) {
      stop("circular variables need to be a vector of variable names or
           indices.")
    }
    toclust1 <- toclust[cir.var]
    if (!is.data.frame(toclust1)) {
      stop("undefined variables selected.")
    }
  }

  ## Switch so we only give off one warning.
  assign(".MonoClustwarn", 0, envir = .GlobalEnv)
  ## Right now, each observation has equal weight. This could be made into an
  ## option.
  weights <- rep(1, nrow(toclust))

  ## Categorical Variable Ordering
  quali_ordered <- NULL

  ## REMOVE: Tan, 4/6/15. These codes haven't show usefulness, but R v1.3.1
  ## changed output of lapply to list, fail running.
  ## ADD: Tan, 9/6/20. Change to dplyr.
  ## See which variables are character and convert them to factors.
  toclust <- purrr::modify_if(toclust, is.character, as.factor)

  ## If we have categorical variables we want to use Gower's distance unless
  ## otherwise specified.
  if (any(purrr::map_lgl(toclust, is.factor)) && is.null(distmethod)) {
    distmethod <- "gower"
  }

  ## Tan 4/16/17, if there is a circular variable in the data set, use Gower's
  ## distance unless otherwise specified.
  if (!is.null(cir.var) && is.null(distmethod)) {
    distmethod <- "gower"
  }

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

  toclust0 <- toclust

  ## Determine whether we need to use PCAmix, or whether we only have
  ## quantitative variables.
  factors <- purrr::map_lgl(toclust, is.factor)
  qualtog <- any(factors)
  quanttog <- any(!factors)

  # If has at least one categorical variable
  extracols <- NULL
  catnames <- character(0)
  # TODO: Candidate to move to a function
  # Input: toclust, factors, corders
  # Output: toclust, quantis
  if (qualtog) {
    ## Since we are going to have multiple factor orderings,
    ## we need to consider each of these orderings a separate "explanatory"
    ## variable, and determine the number of extra columns we will need based on
    ## the number of orderings and the number of factors in the dataset.
    num_factors <- sum(factors)
    extracols <- ncol(toclust) + 1:(num_factors * (corders - 1))

    ## Use PCAmix with quantitative variables the non-factors and qualitative
    ## variables the factors.
    pca <- PCAmixdata::PCAmix(X.quanti = toclust[, !factors],
                              X.quali = toclust[, factors],
                              graph = FALSE)

    ## Get the number of levels of each factor.
    numbyvar <- purrr::map_int(toclust[, factors], nlevels)

    ## Get each factor variable name.
    catnames <- colnames(toclust)[factors]

    ## Sort the levels alphabetically as this is the way they will come out of
    ## PCAmix.
    catrepvarlevel <- purrr::flatten_chr(purrr::map(toclust[, factors],
                                                    ~ sort(levels(.x))))
    names(catrepvarlevel) <- rep(catnames, numbyvar)

    ## Set up variables to be filled in the loop
    ## This method seems a bit odd, but PCAmix looks at all levels of all
    ## variables simultaneously, and some work must be done to keep the levels
    ## with the factor they belong to.
    catrepvarlevelordered <- character(0)
    cuts_quali <- numeric(0)
    position <- 0

    for (j in 1:corders) {
      ## Do this for each dimension that we want 1-5 (number of orderings)
      ## Order and put them in columns
      catrepvarlevel <- catrepvarlevel[order(pca$categ.coord[, j])]
      catrepvarlevelordered <- cbind(catrepvarlevelordered, catrepvarlevel)

      for (p in seq_len(num_factors)) {
        position <- position + 1
        quali_ordered[[position]] <-
          catrepvarlevel[names(catrepvarlevel) == catnames[p]]

        ## We want the orders and the cuts
        ## The cuts are the ordered levels that we can use to partition the
        ## dataset. They are numeric, but correspond to a division between two
        ## levels of a categorical variable according to a particular ordering
        ## from PCAmix.
        fcolumn <- toclust[, which(factors)[p]]
        cuts_quali <-
          cbind(cuts_quali,
                purrr::map_dbl(fcolumn,
                               ~ match(.x, quali_ordered[[position]]) + 1))
      }
    }

    qual_quant <- cuts_quali - 1
    toclust[, factors] <- qual_quant[, 1:num_factors]
    toclust[, extracols] <- qual_quant[, -c(1:num_factors)]

  }
  ## Now, we have extra columns, we have to make sure that we do not consider
  ## them quantitative.
  quantis <- !factors
  quantis <- c(quantis, rep(FALSE, length(extracols)))

  # CLUSTERING ON CIRCULAR VARIABLE
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

  # Clustering on common quantitative variables
  cuts <- toclust
  if (qualtog) cuts[, c(which(factors), extracols)] <- cuts_quali
  if (quanttog) {
    cuts_quant <- purrr::map_dfc(toclust %>%
                                   dplyr::select(which(quantis)),
                                 find_closest)
    cuts[, which(quantis)] <- cuts_quant
  }
  # if (!qualtog) catnames <- character(0)

  ## Variables that are simple derivatives of inputs that will be used a lot.
  labs <- labels
  nobs <- dim(toclust)[1]
  # nvars <- dim(toclust)[2]

  distmats <- matrix()
  # data <- as.data.frame(toclust)

  ## which column to use in the distance matrix
  # distcols <- c(1:ncol(toclust0),rep(which(factors),corders-1))

  ## Since we have multiple categorical orderings for each factor variable,
  ## We need to label each of these orderings. I chose to do it in a kind of
  ## bizarre way.
  ## The name followed by *~* and then the number.
  ## This is odd so that we can replace the *~* in the end which is a unique
  ## character combination and likely wont interfere with any of the variable
  ## names.
  if (length(catnames) > 0) {
    othercolnames <- paste(colnames(toclust)[rep(which(factors), corders - 1)],
                           rep(c(2:corders), each = sum(factors)),
                           sep="*~*")
    currcolnames  <- paste(colnames(toclust)[which(factors)], 1, sep="*~*")
    colnames(toclust)[which(factors)] <- currcolnames
    colnames(toclust)[-c(1:ncol(toclust0))] <- othercolnames
    catnames <- c(currcolnames, othercolnames)
  }

  colnames(cuts) <- colnames(toclust)

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

  distmats <- as.matrix(distmat0)

  members <- seq_len(nobs)

  ## Set up a vector containing each observation's membership. Put into global
  ## environment, but this will be deleted at the end of this function. Using
  ## global environment allows us to modify things recursively as we partition
  ## clusters.
  assign(".Cloc", rep(1,nobs), envir = .GlobalEnv)

  ## Likewise, set up the first (entire dataset) cluster in our Cluster frame
  ## where we keep track of each of the clusters and the partitioning.
  assign(".Cluster_frame",
         data.frame(number = 1,
                    var = "<leaf>",
                    n = nobs,
                    wt = sum(weights[members]),
                    inertia = inertia_calc(distmats[members, members]),
                    bipartvar = "NA",
                    bipartsplitrow = NA,
                    bipartsplitcol = NA,
                    inertiadel = 0,
                    yval = 1,
                    medoid = medoid(members, distmats),
                    category = NA,
                    cut = NA,
                    loc = 0.1,
                    stringsAsFactors = FALSE,
                    split.order = 0),
         envir = .GlobalEnv)

  split.order <- 1
  ## This loop runs until we have nclusters, have exhausted our observations or
  ## run into our minbucket/minsplit restrictions.
  while (sum(.Cluster_frame$var=="<leaf>") < nclusters) {
    check <- checkem(toclust, cuts, distmats, catnames, variables, weights,
                     minsplit, minbucket, split.order)
    split.order <- split.order + 1
    if (check == 0) break
  }

  ## Most of the rest of the function does some bizarre text operations
  ## the reason for this is because I stole a lot of code from rpart,
  ## so we need to follow their text and labelling conventions to that our objects
  ## which inherit from rpart can print and plot correctly.

  ## Change the number column to rownames...
  rownames(.Cluster_frame)<-.Cluster_frame$number
  .Cluster_frame2<-.Cluster_frame[,-1]

  ## This is what will print at each terminal node on the dendrogram
  ## (See plot.MonoClust).
  textfxn<-function(yval,dev,wt,ylevel,digits,n,meds,names, use.n){
    paste("\n  n=", n,"\n  M=",meds, sep="")
  }

  ## Seperate categorical and quantitative splits as the text and plot
  ## functions must treat them a bit differently
  var <-.Cluster_frame2$var
  cattog <- .Cluster_frame2$category

  splits<-which(var != '<leaf>')
  cat_splits<-which(var != '<leaf>' & cattog == 1)

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

  ## MODIFY: Tan, 12/14. Change input of getlevels. If getlevels doesn't see the whole structure of output, we can't
  ## set correct left and right node labels. See new getlevels function for more details.
  # labs<-c('root',sapply(splits,getlevels,cats = cat_splits,varnames=var, frame=.Cluster_frame2,catnames=catnames,quali_ordered=quali_ordered))
  labs<-getlevels(splits, cats = cat_splits,varnames=var, frame=.Cluster_frame2,catnames=catnames,quali_ordered=quali_ordered, digits=digits)

  ## name a column what I probably should hav already named it, but I don't want to change all the code.
  colnames(.Cluster_frame2)[4] <-"dev"

  ## Reorder the columns so they print out nicely, again because I don't want to go back and change things.
  .Cluster_frame2 <- .Cluster_frame2[,c(1,12,2,3,4,5,6,7,8,9,10,11,13,14)]

  ## This follows somewhat odd rpart conventions.
  dendfxns<-list("text"=textfxn)

  ## Take variables' names out of original data set
  Terms <- colnames(toclust0)

  ## ADD: Tan, 12/15, calculate the mean of each cluster
  centroids <- find.centroid(toclust0, qualtog, quanttog)

  ## ADD: Tan, 4/22/16, add medoids of each cluster
  medoids <- .Cluster_frame2[.Cluster_frame2$var =="<leaf>","medoid"]
  names(medoids) <- rownames(.Cluster_frame2[.Cluster_frame2$var =="<leaf>",])

  # We will return a MonoClust object that also inherits from rpart with all of the neccesary components.
  ## MODIFY: Tan, 12/14, I don't know the difference between labels and labelsnum output
  ## although both of them are used in labels.MonoClust. Maybe for categorical variables?
  rpartobj<-list("frame"=.Cluster_frame2,"labels"=labs,"labelsnum" = labs, "functions"=dendfxns,"qualordered"=quali_ordered, Membership =.Cloc, Dist=distmats, Catnames=catnames,
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

## ADD, Tan, 12/15, function to calculate the mean of each cluster.
## Currently do not work for categorical variables
find.centroid <- function(toclust, qualtog, quanttog) {
  # Don't calculate if there is qualitative variable
  if (qualtog) NA

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

getlevels <- function(ind,cats,varnames,frame,catnames,quali_ordered, digits=getOption('digits')){
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
  # If there exists categorical cutpoint
  if (any(ind %in% cats)) {
    sind <- ind[ind %in% cats == 1]
    for (i in sind) {
      name <- varnames[i]
      qualind <- which(catnames==varnames[i])
      lsplit[which(ind == i)] <- paste(quali_ordered[[qualind]][1:(frame$cut[i]-1)],collapse=" ")
      rsplit[which(ind == i)] <- paste(quali_ordered[[qualind]][-c(1:(frame$cut[i]-1))],collapse=" ")
    }

  }

  node <- as.numeric(row.names(frame))
  parent <- match(node%/%2, node[ind])
  odd <- (as.logical(node%%2))

  labels <- character(nrow(frame))
  labels[odd] <- rsplit[parent[odd]]
  labels[!odd] <- lsplit[parent[!odd]]
  labels[1] <- "root"
  labels
}


splitter<-function(splitrow,data,cuts,dist,catnames,weights, split.order = 0){
  ## This function does the actual act of partitioning, given the row that is to be split "splitrow"

  number <- .Cluster_frame$number[splitrow]
  mems   <- which(.Cloc == number)
  split  <- c(.Cluster_frame$bipartsplitrow[splitrow],.Cluster_frame$bipartsplitcol[splitrow])

  datamems<-data.frame(data[mems,])
  cutsmems<-cuts[mems,]


  memsA <- mems[which(datamems[,split[2]] < cutsmems[split[1],split[2]])]
  memsB <-setdiff(mems,memsA)

  # Tan, 2/17, call the permutation test function to test on the newly created group
  # perm.test condition should be uncommented later when successfully tested
  # if (perm.test) {
  # ptest.result <- permtest(split[2], data, memsA, memsB)
  # ptest.result$aov.tab[1,6]
  # }

  # Tan, 10/3, cutpoint is the middle point between two closest points in two clusters
  datamemsA <- datamems[datamems[,split[2]] < cutsmems[split[1],split[2]], split[2]]
  datamemsB <- setdiff(datamems[,split[2]], datamemsA)
  mid.cutpoint <- (max(datamemsA) + min(datamemsB))/2

  ## Make the new clusters.
  Anum<-number*2
  Bnum<-number*2+1

  .Cloc[memsA]<<-Anum
  .Cloc[memsB]<<-Bnum

  variable <- colnames(data)[.Cluster_frame$bipartsplitcol[splitrow]]

  ## This seperates the categorical variable from the level.
  ## Probably bad coding, parsing strings over and over.
  if(grepl(variable,"*~*",fixed=TRUE)){
    variable <- strsplit(variable,"*~*",fixed=TRUE)[[1]]
  }

  ## Is the split categorical?
  if(variable %in% catnames){
    .Cluster_frame[splitrow,12] <<- 1
  }else { .Cluster_frame[splitrow,12] <<- 0 }

  nr<-nrow(.Cluster_frame)

  ## The old cluster now changes some attributes after splitting.
  .Cluster_frame[splitrow,2] <<- variable
  .Cluster_frame[splitrow,6] <<- variable
  #.Cluster_frame[splitrow,13] <<- cutsmems[split[1],split[2]]
  .Cluster_frame[splitrow,13] <<- mid.cutpoint # Use new cutpoint
  .Cluster_frame[splitrow,15] <<- split.order

  ## Tan, 12/14, trying to change the order of rows in Cluster_frame, so that
  ## the children will stay right after the parent, instead of at the very end of the table
  ## Reason: mimic rpart, help print, and plot.
  ## All the comment or new one will be noted carefully so that it can be rollbacked easily

  # ADD: Tan, 12/14.
  # If current row is not the last row in the table, shift all the rows two rows down to make space for the children
  # Otherwise, proceed normally
  if (splitrow < nrow(.Cluster_frame)) {
    emptyrows <- rbind(.Cluster_frame[1,], .Cluster_frame[1,])
    emptyrows[,] <- NA
    .Cluster_frame <<- rbind(head(.Cluster_frame, splitrow),
                             emptyrows,
                             tail(.Cluster_frame, -splitrow))
  }

  ## This part is the same as the commented code, just replace nr by splitrow
  ## New cluster 1 gets some new attributes
  .Cluster_frame[splitrow+1,1] <<- Anum
  .Cluster_frame[splitrow+1,2] <<- "<leaf>"
  .Cluster_frame[splitrow+1,3] <<- length(memsA)
  .Cluster_frame[splitrow+1,4] <<- sum(weights[memsA])
  .Cluster_frame[splitrow+1,5] <<- inertia_calc(dist[memsA,memsA])
  .Cluster_frame[splitrow+1,10] <<- 1-.Cluster_frame[splitrow,9]/.Cluster_frame[1,5]
  .Cluster_frame[splitrow+1,11] <<- medoid(memsA,dist)
  .Cluster_frame[splitrow+1,14] <<- .Cluster_frame[splitrow,14] - 1/nr

  ## As does new cluster 2.
  .Cluster_frame[splitrow+2,1] <<- Bnum
  .Cluster_frame[splitrow+2,2] <<- "<leaf>"
  .Cluster_frame[splitrow+2,3] <<- length(memsB)
  .Cluster_frame[splitrow+2,4] <<- sum(weights[memsB])
  .Cluster_frame[splitrow+2,5] <<- inertia_calc(dist[memsB,memsB])
  .Cluster_frame[splitrow+2,10] <<- 1-.Cluster_frame[splitrow,9]/.Cluster_frame[1,5]
  .Cluster_frame[splitrow+2,11] <<- medoid(memsB,dist)
  .Cluster_frame[splitrow+2,14] <<- .Cluster_frame[splitrow,14] + 1/nr


  #     ## REMOVE: Tan, 12/14 New cluster 1 gets some new attributes
  #     .Cluster_frame[nr+1,1] <<- Anum
  #     .Cluster_frame[nr+1,2] <<- "<leaf>"
  #     .Cluster_frame[nr+1,3] <<- length(memsA)
  #     .Cluster_frame[nr+1,4] <<- sum(weights[memsA])
  #     .Cluster_frame[nr+1,5] <<- inertia_calc(dist[memsA,memsA])
  #     .Cluster_frame[nr+1,10] <<- 1-.Cluster_frame[splitrow,9]/.Cluster_frame[1,5]
  #     .Cluster_frame[nr+1,11] <<- medoid(memsA,dist)
  #     .Cluster_frame[nr+1,14] <<- .Cluster_frame[splitrow,14] - 1/nr
  #
  #     ## As does new cluster 2.
  #     .Cluster_frame[nr+2,1] <<- Bnum
  #     .Cluster_frame[nr+2,2] <<- "<leaf>"
  #     .Cluster_frame[nr+2,3] <<- length(memsB)
  #     .Cluster_frame[nr+2,4] <<- sum(weights[memsB])
  #     .Cluster_frame[nr+2,5] <<- inertia_calc(dist[memsB,memsB])
  #     .Cluster_frame[nr+2,10] <<- 1-.Cluster_frame[splitrow,9]/.Cluster_frame[1,5]
  #     .Cluster_frame[nr+2,11] <<- medoid(memsB,dist)
  #     .Cluster_frame[nr+2,14] <<- .Cluster_frame[splitrow,14] + 1/nr

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
#' @param row The row in .Cluster_frame that would be assessed.
#' @inheritParams checkem
#'
#' @return This function changes the .Cluster_frame in global environment, it's
#'   not supposed to return anything.
#' @importFrom foreach `%dopar%`
#' @keywords internal
find_split <- function(row, data, cuts, dist, variables, weights, minsplit,
                       minbucket) {

  frame <- .Cluster_frame
  bycol <- numeric()
  number <- frame[row, 1]
  mems <- which(.Cloc == number)
  inertiap <- frame[row, 5]

  if (inertiap == 0 | frame[row, 3] < minsplit | frame[row, 3] == 1) {
    # Tan 9/24 This is one obs cluster. Set bipartsplitrow value 0 to stop
    # checkem forever.
    # MG, 9/25 I think this means we won't explore this node again. But make it
    # doesn't stop the search into other nodes.
    # Yup, but we waste resources by keep checking them again and again. The
    # "candidates" of checkem keeps getting longer (at most n) instead of just 2
    # new splits.
    frame[row, 7] <- 0
    .Cluster_frame <<- frame
    return(0)
  }

  ## Subset the data and cut matricies
  # MODIFY: Tan, 7/3/16, add search space limit
  # datamems<-data[mems,]
  # cutsmems<-cuts[mems,]
  datamems <- data.frame(data[mems, variables])
  cutsmems <- data.frame(cuts[mems, variables])

  ## For each possible cut, calculate the inertia. This is where using a
  ## discrete optimization algorithm would help a lot.
  bycol <- foreach::foreach(i = seq_len(ncol(datamems)),
                   .combine = cbind,
                   .exports = c("datamems", "cutsmems", "dist", "mems")) %dopar%
    function(index) {
      data_col <- datamems[, index]
      cuts_col <- cutsmems[, index]

      sapply(cuts_col, function(x) {
        memsA <- mems[which(data_col < x)]
        memsB <- setdiff(mems, memsA)
        ifelse(length(memsA) * length(memsB) == 0,
               NA,
               inertia_calc(dist[memsA, memsA]) +
                 inertia_calc(dist[memsB, memsB]))
      })
    }

  # for (i in seq_len(ncol(datamems))) {
  #
  #   data_col <- datamems[, i]
  #   cuts_col <- cutsmems[, i]
  #
  #   bycol <- cbind(bycol, sapply(cuts_col, function(x) {
  #     memsA <- mems[which(data_col < x)]
  #     memsB <- setdiff(mems, memsA)
  #     ifelse(length(memsA) * length(memsB) == 0,
  #            NA,
  #            inertia_calc(dist[memsA, memsA]) + inertia_calc(dist[memsB, memsB]))
  #     }))
  # }
  # Difference between current cluster and the possible splits
  vals <- inertiap - bycol
  ## Say no diference if we have NA or infinite (happens when no split is possible)
  vals[!is.finite(vals) | is.na(vals)] <- 0

  ## This is the best split.
  maxval <- max(vals)

  ## This is the maximum inertia change indep
  ind <- which((inertiap - bycol) == maxval, arr.ind= TRUE)

  ## Add one more column to check minbucket
  ind <- cbind(ind, TRUE)

  for (i in 1:nrow(ind)) {
    split <- ind[i, ]
    left.size <- sum(datamems[, split[2]] < cutsmems[split[1], split[2]])
    right.size <- length(mems) - left.size
    if (left.size < minbucket | right.size < minbucket) ind[i, 3] <- FALSE
  }

  ## Remove all row doesn't satisfy minbucket and make sure output is always a matrix even though it has only one row
  ind <- matrix(ind[!ind[,3] == FALSE, ], ncol = 3)

  ## If multiple splits produce the same inertia change output a warning.
  #if(nrow(ind) > 1 & .MonoClustwarn==0){.MonoClustwarn <<- 1; warning("One or more of the splits chosen had an alternative split that reduced deviance by the same amount.")}

  # If there is some row that satisfies minbucket
  if (nrow(ind) != 0) {
    split <- ind[1, ]

    memsA <- mems[which(datamems[,split[2]] < cutsmems[split[1], split[2]])]
    memsB <- setdiff(mems, memsA)

    # calculate our change in inertia
    inertiadel <- inertiap -
      inertia_calc(dist[memsA,memsA]) -
      inertia_calc(dist[memsB,memsB])

    ## Update our frame
    frame[row,7] <- split[1]
    frame[row,8] <- variables[split[2]]
    frame[row,9] <- inertiadel
  } else  # Otherwise, stop as a leaf
    frame[row,7] <- 0


  .Cluster_frame <<- frame
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
#' @param dist Distance matrix of all observations in the data.
#' @param catnames Categorical variables' names.
#' @param weights Weights on each observation. Hasn't been implemented in
#' exported function yet. Vector of 1 for all observations.
#' @param minsplit The minimum number of observations that must exist in a node
#'   in order for a split to be attempted.
#' @param split.order The control argument to see how many split has been done.
#' @inheritParams MonoClust
#'
#' @return It is not supposed to return anything because global environment was
#'   used. However, if there is nothing left to split, it returns 0 to tell the
#'   caller to stop running the loop.
#' @keywords internal
checkem <- function(data, cuts, dist, catnames, variables, weights, minsplit,
                    minbucket, split.order = 0) {

  ## Current terminal nodes
  candidates <- which(.Cluster_frame$var == "<leaf>" &&
                        is.na(.Cluster_frame$bipartsplitrow))
  ## Split the best one. Return to nada which never gets output.
  nada <- purrr::map(candidates,
                     ~ find_split(.x, data, cuts, dist, variables, weights,
                                  minsplit, minbucket))
  ## See which ones are left.
  candidates2 <- which(.Cluster_frame$var == "<leaf>" &&
                         .Cluster_frame$bipartsplitrow != 0)
  ## If nothing's left, stop running.
  if (length(candidates2) == 0) return(0)

  ## Find the best inertia change of all that are possible
  maxone <- max(.Cluster_frame$inertiadel[candidates2], na.rm = TRUE)
  splitrow <- candidates2[which(.Cluster_frame$inertiadel[candidates2] ==
                                  maxone)]

  ## Make new clusters from that cluster
  # splitter(splitrow, data, cuts,dist,catnames,weights)
  splitter(splitrow[1], data, cuts, dist, catnames, weights, split.order)
  # Tan, 9/24, in case there are more than one node equal to max
  # MG, 9/25, I thought the earlier code would make sure only one is identified
  # as top but that might not be true. It never caused a problem before.
  # Maybe because we fixed inertia fn, equal inertias occurred for small
  # clusters
}
