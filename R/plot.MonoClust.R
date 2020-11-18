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
#'   average branch length. This parameter is ignored if `uniform = TRUE`.
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
#' @param col.type When `cols` is set, choose whether the color indicators are
#'   shown in a form of solid lines below the leaves (`"l"`), or big points
#'   (`"p"`), or both (`"b"`).
#' @param rel.loc.x Whether to use the relative distance between clusters as x
#'   coordinate of the leaves. Default is TRUE.
#' @param show.pval If MonClust object has been run through [perm.test()],
#'   whether to show p-value on the tree.
#' @param ... Arguments to be passed to [graphics::plot.default()] and
#'   [graphics::lines()].
#' @inheritParams print.MonoClust
#'
#' @return A plot of splitting rule.
#' @export
#'
#' @examples
#' library(cluster)
#' data(ruspini)
#'
#' # MonoClust tree
#' ruspini4sol <- MonoClust(ruspini, nclusters = 4)
#' plot(ruspini4sol)
#' \donttest{
#' # MonoClust tree after permutation test is run
#' ruspini6sol <- MonoClust(ruspini, nclusters = 6)
#' ruspini6_test <- perm.test(ruspini6sol,
#'                            data = ruspini,
#'                            method = "sw",
#'                            rep = 1000)
#' plot(ruspini6_test, branch = 1, uniform = TRUE)
#' }
plot.MonoClust <- function(x, uniform = FALSE, branch = 1,
                           margin = c(0.12, 0.02, 0, 0.05),
                           minbranch = 0.3, text = TRUE, which = 4,
                           stats = TRUE,
                           abbrev = c("no", "short", "abbreviate"),
                           digits = getOption("digits") - 2,
                           cols = NULL, col.type = c("l", "p", "b"),
                           rel.loc.x = TRUE, show.pval = TRUE, ...) {

  if (!inherits(x, "MonoClust"))
    stop("Not a MonoClust object")
  if (!(which %in% 1:4))
    stop("\"which\" has to be a value between 1 and 4.")
  abbrev <- match.arg(abbrev)
  if (!is.null(cols)) {
    if (length(cols) > 1) {
      if (length(cols) < sum(x$frame$var == "<leaf>"))
        stop("When set, \"col\" has to contain 1 color or the number of colors
        is greater or equal to the number of leaves.")
      cols <- cols[1:sum(x$frame$var == "<leaf>")]
    }

    col.type <- match.arg(col.type)
  }

  plot_tree(x, uniform = uniform, branch = branch, margin = margin,
            minbranch = minbranch, rel.loc.x = rel.loc.x, ...)

  if (text) {
    text_tree(x, which = which, digits = digits, stats = stats, abbrev = abbrev,
              cols = cols, cols.type = col.type, rel.loc.x = rel.loc.x,
              show.pval = show.pval, uniform = uniform, minbranch = minbranch)

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
