test_that("ggcv checks if the object is a well-formed cv.MonoClust object", {

  example_tab <-
    # Wrong variable name (Std.Dev.)
    tibble::tibble(ncluster = 2:4,
                   MSE = c(9721, 3843, 3923),
                   `Std.Dev.` = c(2570, 1000, 1000))

  expect_error(
    ggcv(),
    "\"cv.obj\" is required."
  )

  expect_error(
    ggcv(example_tab),
    "Not a legitimate \"cv.MonoClust\" object.")

  example_cv <- list(cv = example_tab)
  class(example_cv) <- "cv.MonoClust"

  expect_error(
    ggcv(example_cv),
    "Not a legitimate \"cv.MonoClust\" object.")
})

test_that("ggcv returns correct plot", {
  example_tab <-
    tibble::tibble(ncluster = 2:4,
                   MSE = c(9721, 3843, 3923),
                   `Std. Dev.` = c(2570, 1000, 1000))
  example_cv <- list(cv = example_tab)
  class(example_cv) <- "cv.MonoClust"
  p <- ggcv(example_cv)

  expect_s3_class(p$layers[[1]], "ggproto")
  expect_s3_class(p$layers[[2]]$geom, "GeomPoint")
  expect_s3_class(p$layers[[3]]$geom, "GeomLine")

  p1 <- ggcv(example_cv, type = "b")
  expect_identical(p1, p)

  p2 <- ggcv(example_cv, type = "l")
  expect_s3_class(p2$layers[[1]], "ggproto")
  expect_s3_class(p2$layers[[2]]$geom, "GeomLine")
  expect_error(p2$layers[[3]]$geom)

  p3 <- ggcv(example_cv, type = "p")
  expect_s3_class(p3$layers[[1]], "ggproto")
  expect_s3_class(p3$layers[[2]]$geom, "GeomPoint")
  expect_error(p3$layers[[3]]$geom)

})
