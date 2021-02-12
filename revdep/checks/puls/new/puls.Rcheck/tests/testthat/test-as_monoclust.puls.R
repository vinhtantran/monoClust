test_that("print error when coercing wrong PULS class to MonoClust class", {
  testthat::expect_warning({
    a <- list(test = 1:3)
    as_MonoClust(a)
  },
  "as_MonoClust does not know how to handle object of class  list .")

  testthat::expect_error({
    b <- list(test = 1:3, y = 5)
    class(b) <- "PULS"
    as_MonoClust(b)
  },
  "Object needs a \"frame\". See PULS.object for details.")

  testthat::expect_error({
    b <- list(frame = "ABC",
              test = 1:3, y = 5, membership = rep(1, 10),
              dist = matrix(5, 2, 2),
              terms = c("a", "b"),
              medoids = c(10, 15))
    class(b) <- "PULS"
    as_MonoClust(b)
  },
  "Object must have required objects. See PULS.object for details.")

  testthat::expect_error({
    a <- list(test = 1:3, frame = "ABC",
              test = 1:3, y = 5, membership = rep(1, 10),
              dist = matrix(5, 2, 2),
              terms = c("a", "b"),
              medoids = c(10, 15),
              alt = FALSE)
    class(a) <- "PULS"
    as_MonoClust(a)
  },
  "\"frame\" object must be a data.frame or a data.frame derivation.")

  testthat::expect_error({
    a <- list(test = 1:3, frame = data.frame(var1 = 1, var2 = "string"),
              y = 5, membership = rep(1, 10),
              dist = matrix(5, 2, 2),
              terms = c("a", "b"),
              medoids = c(10, 15),
              alt = FALSE)
    class(a) <- "PULS"
    as_MonoClust(a)
  },
  "\"frame\" must have required columns. See PULS.object for details.")

})

test_that("correctly coerce PULS class to MonoClust class", {
  testthat::expect_s3_class({
    a <- list(frame = new_node(number = 1,
                               var = "<leaf>",
                               n = 30,
                               wt = 10,
                               bipartsplitrow = 5,
                               bipartsplitcol = 10,
                               inertia = 10000,
                               inertia_explained = 0.5,
                               medoid = 2,
                               loc = 0.5,
                               alt = FALSE),
              membership = rep(1, 10),
              dist = matrix(5, 2, 2),
              terms = c("a", "b"),
              medoids = c(10, 15),
              alt = FALSE)
    class(a) <- "PULS"
    as_MonoClust(a)
  },
  "MonoClust")
})
