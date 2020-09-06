
labels.MonoClust <- function(object, abbrev, ...) {
    if (abbrev == 0) {
        labs <- gsub("*~*", ".", object$labelsnum, fixed = TRUE)
        return(labs)
    } else if (abbrev == "L") {
        vars <- object$frame$var
        uvars <- unique(vars)
        names <- uvars[uvars != "<leaf>"]
        nums <- paste("V", 1:length(names), sep = "")
        
        labs <- object$labelsnum
        subit <- function(name, num, VECT) sapply(VECT, function(x) gsub(name, num, x, fixed = TRUE))
        apply(cbind(names, nums), 1, function(row) labs <<- subit(row[1], row[2], labs))
        return(labs)
        
        
    } else {
        getabbrevs <- sapply(object$labels, abbreviate.t, abbrev = abbrev)
        return(getabbrevs)
    }
}








