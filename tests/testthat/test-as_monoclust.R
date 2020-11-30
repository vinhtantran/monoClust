test_that("as_MonoClust shows correct warning when called", {
  testthat::expect_warning(as_MonoClust(iris),
               "as_MonoClust does not know how to handle object of class")
})
