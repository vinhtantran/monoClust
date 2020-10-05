#' Print Monothetic Clustering Results
#'
#' Render the `MonoClust` split tree in an easy to read format with important
#' information such as terminal nodes, p-value (if possible), etc.
#'
#' @param x MonoClust result object.
#' @param abbrev Whether to print the abbreviated versions of variable names.
#'   Can be either "no" (default), "short", or "abbreviate". Short forms of them
#'   can also be used.
#'
#'   If "no", the labels recorded in `x$labels` are used.
#'
#'   If "short", variable names will be turned into "V1", "V2", ...
#'
#'   If "abbreviate", [abbreviate()] function will be used. Use the optional
#'   arguments for this function.
#' @param spaces Spaces indent between 2 tree levels.
#' @param digits Number of significant digits to print.
#' @param ... Optional arguments to [abbreviate()]
#'
#' @return A nicely displayed MonoClust split tree.
#' @export
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#' ruspini4sol <- MonoClust(ruspini, nclusters = 4)
#' print(ruspini4sol)
print.MonoClust <- function(x, abbrev = c("no", "short", "abbreviate"),
                            spaces = 2L, digits = getOption("digits"), ...) {

  if (!inherits(x, "MonoClust"))
    stop("Not a legitimate \"MonoClust\" object")

  abbrev <- match.arg(abbrev)

  frame <- x$frame
  # REMOVE: Tan, 9/14/20. Not necessary
  # ylevel <- attr(x, "ylevels")
  node <- frame$number
  depth <- tree_depth(node)
  # indent <- stringr::str_pad("", 32L)

  # 2L is because of 1 number (1 or 2 digits) and the bracket
  indent <- stringr::str_pad(stringr::str_c(node, ")"), depth * spaces + 2L)

  # MODIFY: Tan, 9/14/20. Not necessary
  # tfun <- (x$functions)$print
  # if (!is.null(tfun)) {
  #   if (is.null(frame$yval2))
  #     yval <- tfun(frame$yval, ylevel, digits) else yval <- tfun(frame$yval2,
  # ylevel, digits)
  # } else yval <- format(signif(frame$yval, digits = digits))

  yval <- format(signif(frame$yval, digits = digits))

  # REMOVE: Tan, 12/14. This line seems to be a debug line. Mo meaning.
  ## print(yval)
  ## ADD: Tan, 3/1/15. Add p value check.
  has_pvalue <- !is.null(frame$p.value)
  term <- rep(" ", length(depth))
  term[frame$var == "<leaf>"] <- "*"
  labs <- create_labels(x, abbrev = abbrev, ...)
  n <- frame$n
  ## MODIFY: Tan, 3/1/15. Add p value.
  ## z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)),
  ## round((1-as.numeric(yval)/1),digits=2), term)

  z <- paste(indent, labs, n,
             format(signif(frame$inertia, digits = digits)),
             format(signif((1 - as.numeric(yval) / 1), digits = digits)),
             ifelse(has_pvalue, frame$p.value, ""), term)

  # MODIDY: Tan, 9/14/20. Not necessary.
  # omit <- x$na.action
  # if (length(omit)) {
  #   cat("n=", n[1L], " (", naprint(omit), ")\n\n", sep = "")
  # } else {
  #   cat("n=", n[1L], "\n\n")
  # }

  cat("n=", n[1L], "\n\n")

  cat("Node) Split, N, Cluster Inertia, Proportion Inertia Explained,",
      ifelse(has_pvalue, "Bonferroni adj. p-value", ""), "\n")
  cat("      * denotes terminal node\n\n")
  cat(z, sep = "\n")
  if (!is.null(x$circularroot$var)) {
    cat("Circular variable's first cut\n")
    # cat("Circular variable(s)' first cut\n")
    for (i in seq_len(length(x$circularroot$var))) {
      cat(x$terms[x$circularroot$var[i]], ": ", x$circularroot$cut[i], "\n")
    }
  }
  return(invisible(x))
}

#' Create Labels for Split Variables
#'
#' This function prints variable's labels for a `MonoClust` tree.
#'
#' @param object MonoClust result object.
#' @inheritParams print.MonoClust
#'
#' @return A named vector of labels corresponding to variable's names (at vector
#'   names).
#' @keywords internal
create_labels <- function(object, abbrev, ...) {
  if (abbrev == "no") {
    # labs <- gsub("*~*", ".", object$labelsnum, fixed = TRUE)
    labs <- object$labels
  } else if (abbrev == "short") {
    vars <- object$frame$var
    uvars <- unique(vars)
    names <- uvars[uvars != "<leaf>"]
    nums <- stringr::str_c("V", seq_len(length(names)))
    labs <- nums
    names(labs) <- names
  } else {
    labs <- purrr::map_chr(object$labels, abbreviate, ...)
  }

  return(labs)
}
