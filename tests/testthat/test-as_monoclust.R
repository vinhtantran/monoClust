test_that("print error when coercing wrong class to MonoClust class", {
  testthat::expect_warning({
    a <- list(test = 1:3)
    as_MonoClust(a)
  },
  "as_MonoClust does not know how to handle object of class list.")
})
