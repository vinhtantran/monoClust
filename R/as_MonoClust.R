#' Coerce Similar Object to MonoClust
#'
#' The function turns a MonoClust-similar object into MonoClust object so it
#' can use supported functions for MonoClust such as [print.MonoClust()] and
#' [plot.MonoClust()].
#'
#' `as_MonoClust()` is an S3 generic. The function itself doesn't run unless
#' it is implemented for another similar object. Currently, this function is not
#' implemented within `monoClust` package.
#'
#' @param x An object that can be coerced to MonoClust object.
#' @param ... For extensibility.
#'
#' @export
as_MonoClust <- function(x, ...) {
  UseMethod("as_MonoClust")
}

#' @export
#' @rdname as_MonoClust
as_MonoClust.default <- function(x, ...) {
  warning(paste("as_MonoClust does not know how to handle object of class ",
                class(x), "."))
}
