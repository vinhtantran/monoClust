#' Plot MonoClust splitting rule tree
#'
#' @param x MonoClust object
#' @param margin Margins to the border of the graphics machine
#' @param which Fill in later
#' @param abbrev Whether abbreviations of labels are used
#' @param text Fill in later
#' @param cols Whether to use color or not
#' @param ... Other args
#'
#' @return Plot of splitting rule
#' @export
#'
#' @examples NULL
plot.MonoClust<-function(x,margin,which,abbrev=4,text=TRUE,cols=NULL,...){
    ## This function sets some defaults and changes things a bit, but is mostly a
    ## wrapper for our slightly modified version of rpart's plot function (see plots.R).

    if(missing(margin)){margin<-c(.12,.02,0,.05)}
    if(missing(which)){which <- 4}

    plot.rpart(x,margin=margin,...)

    ## REMOVE: Tan, 3/1/15, Remove Inertia line
#     lines(x=c(.88,.88),y=c(0,1))

#     for(i in seq(0,1,.1)){
#         lines(x=c(.86,.88),y=c(i,i))
#         text(.73,i,i)
#     }

    if(text){
        text.MonoClust(x,which=which,abbrev=abbrev,cols = cols)
    }

    if (!is.null(x$circularroot$var)) {
      text(x=1, y=1, "Circ root")
      for (i in 1:length(x$circularroot$var)) {
        text(x=1, y=1-i*0.05, paste(x$terms[x$circularroot$var[i]], ": ", x$circularroot$cut[i]))
      }
    }
}

Nclustplot<-function(x, main, type, ylab, xlab,...){

    if(missing(main)){main<-"Marginal Cluster Analysis"}
    if(missing(type)){type<-"b"}
    if(missing(ylab)){ylab<-"Proportion of Deviance Explained"}
    if(missing(xlab)){xlab<-"Number of Clusters"}

    inds <- seq(from=2,to=nrow(x$frame), by =2)
    plot(inds,round((1-as.numeric(x$frame$yval[inds])/1),digits=2),type=type, xaxt="n", ylab=ylab, xlab=xlab, main=main)
    axis(1, at=inds, labels= as.character(2 + 0:(length(inds)-1)))

}



