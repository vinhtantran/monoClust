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
#' @seealso [abbreviate()]
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
  labs <- create_labels(x, abbrev = abbrev, digits = digits, ...)$labels
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

  if (any(frame$alt))
    cat("\nNote: One or more of the splits chosen had an alternative split that
        reduced inertia by the same amount. See \"alt\" column of \"frame\"
        object for details.")

  return(invisible(x))
}

#' Create Labels for Split Variables
#'
#' This function prints variable's labels for a `MonoClust` tree.
#'
#' @inheritParams print.MonoClust
#'
#' @return A list containing two elements:
#'   * `varnames`: A named vector of labels corresponding to variable's names
#'   (at vector names).
#'   * `labels`: Vector of labels of splitting rules to be displayed.
#' @seealso [abbreviate()]
#' @keywords internal
create_labels <- function(x, abbrev, digits = getOption('digits'), ...) {

  frame <- x$frame

  # Rename variable as indicated in abbrev
  vars <- frame$var
  uvars <- unique(vars)
  names <- uvars[uvars != "<leaf>"]

  if (abbrev == "short") {
    varnames <- stringr::str_c("V", seq_len(length(names)))
  } else if (abbrev == "abbreviate") {
    varnames <- purrr::map_chr(names, abbreviate, ...)
  } else {
    varnames <- names
  }

  names(varnames) <- names

  # Create split labels
  split_index <- which(frame$var != "<leaf>")
  lsplit <- rsplit <- character(length(split_index))
  # If there exists quantitative cutpoint
  # if (any(!ind %in% cats)) {
  # sind <- ind[ind %in% cats == 0]
  label <- varnames[frame$var[split_index]]
  level <- frame$cut[split_index]

  # lsplit[ind %in% cats == 0] <- paste(name,"<",round(level, digits),sep=" ")
  # rsplit[ind %in% cats == 0] <- paste(name,">=",round(level, digits),sep=" ")
  lsplit <- paste(label, "<", round(level, digits), sep=" ")
  rsplit <- paste(label, ">=", round(level, digits), sep=" ")
  # }
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

  node <- frame$number
  parent <- match(node %/% 2, node[split_index])
  odd <- as.logical(node %% 2)

  labels <- character(nrow(frame))
  labels[odd] <- rsplit[parent[odd]]
  labels[!odd] <- lsplit[parent[!odd]]
  labels[1] <- "root"

  return(list(varnames = varnames, labels = labels))
}
