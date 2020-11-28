#' Test If The Object is A MonoClust
#'
#' This function returns `TRUE` for MonoClust, and FALSE for all other objects.
#'
#' @param mono_obj An object.
#'
#' @return `TRUE` if the object inherits from the `MonoClust` class.
#' @export
is_MonoClust <- function(mono_obj) {
  if (inherits(mono_obj, "MonoClust")) return(TRUE) else return(FALSE)
}
