#' Print MonoClust Cross-Validation Result#'
#'
#' @param x A `cv.MonoClust` object (output of [cv.test()]).
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(cluster)
#' data(ruspini)
#'
#' # 10-fold cross-validation
#' cp_table <- cv.test(ruspini, minnodes = 2, maxnodes = 4)
#' print(cp_table)
#' }
print.cv.MonoClust <- function(x, ...) {

  if (missing(x))
    stop("\"x\" is required.")
  if (!inherits(x, "cv.MonoClust"))
    stop("Not a legitimate \"cv.MonoClust\" object.")

  cat(x[["cv.type"]], "on a MonoClust object \n\n")
  print(x[["cv"]])
}
