## This file is a hack to remove false flags in R CMD check when using foreach
## iterations.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("iter", "k"))
}
