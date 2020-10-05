#' Plot MonoClust Splitting Rule Tree
#'
#' Print the MonoClust tree in the form of dendrogram.
#'
#' @param x MonoClust result object.
#' @param uniform If TRUE, uniform vertical spacing of the nodes is used; this
#'   may be less cluttered when fitting a large plot onto a page. The default is
#'   to use a non-uniform spacing proportional to the inertia in the fit.
#' @param branch Controls the shape of the branches from parent to child node.
#'   Any number from 0 to 1 is allowed. A value of 1 gives square shouldered
#'   branches, a value of 0 give V shaped branches, with other values being
#'   intermediate.
#' @param margin An extra fraction of white space to leave around the borders of
#'   the tree. (Long labels sometimes get cut off by the default computation).
#' @param minbranch Set the minimum length for a branch to `minbranch` times the
#'   average branch length. This parameter is ignored if `uniform=TRUE`.
#'   Sometimes a split will give very little improvement, or even no improvement
#'   at all. A tree with branch lengths strictly proportional to improvement
#'   leaves no room to squeeze in node labels.
#' @param text Whether to print the labels on the tree.
#' @param which Fill in later
#' @param cols Whether to use color on the labels or not. Only works when `text`
#'   is `TRUE`.
#' @param rel.loc.x Whether to use the relative distance between clusters as x
#'   coordinate of the leaves. Default is TRUE.
#' @param ... Arguments to be passed to [graphics::plot.default()] and
#'   [graphics::lines()]
#' @inheritParams print.MonoClust
#'
#' @return A plot of splitting rule.
#' @export
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#' ruspini4sol <- MonoClust(ruspini, nclusters = 4)
#' plot(ruspini4sol)
plot.MonoClust <- function(x, uniform = FALSE, branch = 1, margin = 0,
                           minbranch = 0.3, text = TRUE, which,
                           abbrev = c("no", "short", "abbreviate"),
                           cols = NULL, rel.loc.x = TRUE, ...) {
  ## This function sets some defaults and changes things a bit, but is mostly
  ## a wrapper for our slightly modified version of rpart's plot function (see
  ## plots.R).

  if (!inherits(x, "MonoClust"))
    stop("Not an MonoClust object")

  if (missing(margin)) {
    margin <- c(0.12, 0.02, 0, 0.05)
  }
  if (missing(which)) {
    which <- 4
  }

  plot_tree(x, uniform = uniform, branch = branch, margin = margin,
            minbranch = minbranch, rel.loc.x = rel.loc.x, ...)

  ## REMOVE: Tan, 3/1/15, Remove Inertia line lines(x=c(.88,.88),y=c(0,1))

  # for(i in seq(0,1,.1)) { lines(x=c(.86,.88),y=c(i,i)) text(.73,i,i) }

  if (text) {
    text_tree(x, which = which, abbrev = abbrev, cols = cols,
              rel.loc.x = rel.loc.x)

    if (!is.null(x$circularroot$var)) {
      graphics::text(x = 1, y = 1, "Circ root")
      for (i in seq_len(length(x$circularroot$var))) {
        graphics::text(x = 1, y = 1 - i * 0.05,
             paste(x$terms[x$circularroot$var[i]], ": ",
                   x$circularroot$cut[i]))
      }
    }
  }
}
