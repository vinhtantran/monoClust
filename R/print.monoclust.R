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
#' @param ... Optional arguments to [abbreviate()].
#'
#' @return A nicely displayed MonoClust split tree.
#' @seealso [abbreviate()]
#' @export
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#' ruspini4sol <- MonoClust(ruspini, nclusters = 4)
#' print(ruspini4sol, digits = 2)
print.MonoClust <- function(x, abbrev = c("no", "short", "abbreviate"),
                            spaces = 2L, digits = getOption("digits"), ...) {

  if (!is_MonoClust(x))
    stop("Not a legitimate \"MonoClust\" object.")

  abbrev <- match.arg(abbrev)

  frame <- x$frame
  node <- frame$number
  depth <- tree_depth(node)

  # 2L is because of 1 number (1 or 2 digits) and the bracket
  indent <- stringr::str_pad(stringr::str_c(node, ")"), depth * spaces + 2L)

  inertia_explained <- ifelse(!is.na(frame$inertia_explained),
                              format(signif(frame$inertia_explained,
                                            digits = digits)),
                              "")
  # Add p value check.
  has_pvalue <- !is.null(frame[["p.value"]])
  term <- rep(" ", length(depth))
  term[frame$var == "<leaf>"] <- "*"
  labs <- create_labels(x, abbrev = abbrev, digits = digits, ...)$labels
  n <- frame$n

  z <- paste(indent, labs, n,
             format(signif(frame$inertia, digits = digits)),
             inertia_explained,
             ifelse(has_pvalue, frame$p.value, ""), term)

  cat("n =", n[1L], "\n\n")

  cat("Node) Split, N, Cluster Inertia, Proportion Inertia Explained,",
      ifelse(has_pvalue, "Bonferroni adj. p-value", ""), "\n")
  cat("      * denotes terminal node\n\n")
  cat(z, sep = "\n")
  if (!is.null(x$circularroot$var)) {
    cat("Circular variable's first cut\n")
    for (i in seq_len(length(x$circularroot$var))) {
      cat(x$terms[x$circularroot$var[i]], ": ", x$circularroot$cut[i], "\n")
    }
  }

  # Alternate route notice
  if (x$alt)
    cat(paste("\nNote: One or more of the splits chosen had an alternative",
              "split that reduced inertia by the same amount. See \"alt\"",
              "column of \"frame\" object for details."))

  return(invisible(x))
}
