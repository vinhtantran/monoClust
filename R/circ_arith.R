#' Add/Subtract Circular Values in Degrees/Radian
#'
#' Add/subtract two circular variables in degrees (`%cd+%` and `%cd-%`) and
#' radian (`%cr+%` and `%cr-%`).
#'
#' @param x,y Circular values in degrees/radians.
#'
#' @return A value between [0, 360) in degrees or [0, 2*pi) in radian.
#' @name circ_arith
#' @examples
#' 90 %cd+% 90
#'
#' 250 %cd+% 200
#'
#' 25 %cd-% 80
#'
#' pi %cr+% (pi/2)
#'
NULL

#' @export
#' @rdname circ_arith
`%cd+%` <- function(x, y) {
  return((x + y) %% 360)
}

#' @export
#' @rdname circ_arith
`%cd-%` <- function(x, y) {
  return((x - y) %% 360)
}

#' @export
#' @rdname circ_arith
`%cr+%` <- function(x, y) {
  return((x + y) %% (2 * pi))
}

#' @export
#' @rdname circ_arith
`%cr-%` <- function(x, y) {
  return((x - y) %% (2 * pi))
}
