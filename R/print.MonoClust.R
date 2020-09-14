#' Print Monothetic Clustering Results
#'
#' @param x MonoClust result object.
#' @param abbrev Whether print the abbreviated version of labels.
#' @param spaces Spaces indent between 2 tree levels.
#' @param digits Number of significant digits to print.
#'
#' @return Splitting tree displayed nicely.
#' @export
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#' ruspini4sol <- MonoClust(ruspini, nclusters = 4)
#' print(ruspini4sol)
print.MonoClust <- function(x, abbrev = 0, spaces = 2L,
                            digits = getOption("digits")) {

  ## Modified from print.rpart.

  minlength <- 0
  if (!inherits(x, "MonoClust"))
    stop("Not a legitimate \"MonoClust\" object")

  frame <- x$frame
  # REMOVE: Tan, 9/14/20. Not necessary
  # ylevel <- attr(x, "ylevels")
  node <- frame$number
  depth <- tree_depth(node)
  indent <- stringr::str_pad("", 32L)

  # 2L is because of 1 number (1 or 2 digits) and the bracket
  indent <- stringr::str_pad(stringr::str_c(node, ")"), depth * spaces + 2L)

  # MODIFY: Tan, 9/14/20. Not necessary
  # tfun <- (x$functions)$print
  # if (!is.null(tfun)) {
  #   if (is.null(frame$yval2))
  #     yval <- tfun(frame$yval, ylevel, digits) else yval <- tfun(frame$yval2, ylevel, digits)
  # } else yval <- format(signif(frame$yval, digits = digits))

  yval <- format(signif(frame$yval, digits = digits))

  # REMOVE: Tan, 12/14. This line seems to be a debug line. Mo meaning.
  ## print(yval)
  ## ADD: Tan, 3/1/15. Add p value check.
  has_pvalue <- ifelse(is.null(frame$p.value), FALSE, TRUE)
  term <- rep(" ", length(depth))
  term[frame$var == "<leaf>"] <- "*"
  z <- create_labels(x, abbrev = abbrev)
  n <- frame$n
  ## MODIFY: Tan, 3/1/15. Add p value.
  ## z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)),
  ## round((1-as.numeric(yval)/1),digits=2), term)
  if (has_pvalue) {
    z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)),
               round((1 - as.numeric(yval)/1), digits = digits),
               frame$p.value, term)
  } else {
    z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)),
               round((1 - as.numeric(yval)/1), digits = digits),
               term)
  }

  # MODIDY: Tan, 9/14/20. Not necessary.
  # omit <- x$na.action
  # if (length(omit)) {
  #   cat("n=", n[1L], " (", naprint(omit), ")\n\n", sep = "")
  # } else {
  #   cat("n=", n[1L], "\n\n")
  # }

  cat("n=", n[1L], "\n\n")

  cat("Node, N, Within Cluster Deviance, Proportion Deviance Explained,",
      ifelse(has_pvalue, "Bonferroni adj. p-value", ""), "\n")
  cat("      * denotes terminal node\n\n")
  cat(z, sep = "\n")
  if (!is.null(x$circularroot$var)) {
    cat("Circular variable(s)' first cut\n")
    for (i in seq_len(length(x$circularroot$var))) {
      cat(x$terms[x$circularroot$var[i]], ": ", x$circularroot$cut[i], "\n")
    }
  }
  return(invisible(x))
}

create_labels <- function(object, abbrev) {
  if (abbrev == 0) {
    labs <- gsub("*~*", ".", object$labelsnum, fixed = TRUE)
    return(labs)
  } else if (abbrev == "L") {
    vars <- object$frame$var
    uvars <- unique(vars)
    names <- uvars[uvars != "<leaf>"]
    nums <- paste("V", 1:length(names), sep = "")

    labs <- object$labelsnum
    subit <- function(name, num, VECT)
      sapply(VECT, function(x) gsub(name, num, x, fixed = TRUE))
    apply(cbind(names, nums), 1, function(row)
      labs <<- subit(row[1], row[2], labs))
    return(labs)

  } else {
    getabbrevs <- sapply(object$labels, abbreviate.t, abbrev = abbrev)
    return(getabbrevs)
  }
}
