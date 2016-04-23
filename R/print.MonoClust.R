print.MonoClust <- function (x, abbrev=0,  spaces = 2, digits = options('digits')$digits, 
                             ...) 
{
    
    ## Modified from print.rpart.
    
    minlength <-0
    if (!inherits(x, "rpart")) 
        stop("Not legitimate rpart object")
    if (!is.null(x$frame$splits)) 
        x <- rpconvert(x)
    
    frame <- x$frame
    ylevel <- attr(x, "ylevels")
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32L), collapse = "")
    if (length(node) > 1L) {
        indent <- substring(indent, 1L, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")", 
                        sep = "")
    }
    else indent <- paste(format(node), ")", sep = "")
    tfun <- (x$functions)$print
    if (!is.null(tfun)) {
        if (is.null(frame$yval2)) 
            yval <- tfun(frame$yval, ylevel, digits)
        else yval <- tfun(frame$yval2, ylevel, digits)
    }
    else yval <- format(signif(frame$yval, digits = digits))
    
    ## REMOVE: Tan, 12/14. This line seems to be a debug line. Mo meaning.
    # print(yval)
    ## ADD: Tan, 3/1/15. Add p value check.
    has.pvalue <- ifelse(is.null(frame$p.value), FALSE, TRUE)
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels(x, digits = digits, abbrev= abbrev, minlength = minlength, ...)
    n <- frame$n
    ## MODIFY: Tan, 3/1/15. Add p value.
    # z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)), round((1-as.numeric(yval)/1),digits=2), term)
    if (has.pvalue) {
        z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)), round((1-as.numeric(yval)/1),digits=2), frame$p.value,term)
    } else {
        z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)), round((1-as.numeric(yval)/1),digits=2), term)        
    }
    
    omit <- x$na.action
    if (length(omit)) 
        cat("n=", n[1L], " (", naprint(omit), ")\n\n", sep = "")
    else cat("n=", n[1L], "\n\n")
    cat("Node, N, Within Cluster Deviance, Proportion Deviance Explained,", ifelse(has.pvalue, "Bonferroni adj. p-value", ""),"\n")
    cat("      * denotes terminal node\n\n")
    cat(z, sep = "\n")
    return(invisible(x))
}

