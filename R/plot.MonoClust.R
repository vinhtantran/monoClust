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
#' @param which Labeling modes, which are:
#'   * 1: only splitting variable names are shown, no splitting rules.
#'   * 2: only splitting rules to the left branches are shown.
#'   * 3: only splitting rules to the right branches are shown.
#'   * 4 (default): splitting rules are shown on both sides of branches.
#' @param stats Whether to show statistics (cluster sizes and medoid points) on
#'   the tree.
#' @param cols Whether to shown color bars at leaves or not. It helps matching
#'   this tree plot with other plots whose cluster membership were colored. It
#'   only works when `text` is `TRUE`. Either `NULL`, a vector of one color, or
#'   a vector of colors matching the number of leaves.
#' @param cols_type When `cols` is set, choose whether the color indicators are
#'   shown in a form of solid lines below the leaves (`"l"`), or big points
#'   (`"p"`), or both (`"b"`).
#' @param rel_loc_x Whether to use the relative distance between clusters as x
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
plot.MonoClust <- function(x, uniform = FALSE, branch = 1,
                           margin = c(0.12, 0.02, 0, 0.05),
                           minbranch = 0.3, text = TRUE, which = 4,
                           stats = TRUE,
                           abbrev = c("no", "short", "abbreviate"),
                           digits = getOption("digits") - 2,
                           cols = NULL, cols_type = c("l", "p", "b"),
                           rel_loc_x = TRUE, ...) {
  ## This function sets some defaults and changes things a bit, but is mostly
  ## a wrapper for our slightly modified version of rpart's plot function (see
  ## plots.R).

  if (!inherits(x, "MonoClust"))
    stop("Not a MonoClust object")
  if (!(which %in% 1:4))
    stop("\"which\" has to be a value between 1 and 4.")
  abbrev <- match.arg(abbrev)
  if (!is.null(cols)) {
    if (length(cols) > 1 & length(cols) != sum(x$frame$var == "<leaf>"))
      stop("When set, \"col\" has to contain 1 color or number of colors equal
      to the number of leaves.")
    cols_type <- match.arg(cols_type)
  }

  plot_tree(x, uniform = uniform, branch = branch, margin = margin,
            minbranch = minbranch, rel_loc_x = rel_loc_x, ...)

  ## REMOVE: Tan, 3/1/15, Remove Inertia line lines(x=c(.88,.88),y=c(0,1))

  # for(i in seq(0,1,.1)) { lines(x=c(.86,.88),y=c(i,i)) text(.73,i,i) }

  if (text) {
    text_tree(x, which = which, digits = digits, stats = stats, abbrev = abbrev,
              cols = cols, cols_type = cols_type, rel_loc_x = rel_loc_x)

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
