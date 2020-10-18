## This file is a hack to remove false flags in R CMD check when using foreach
## iterations.
utils::globalVariables(c("iter", "k"))
