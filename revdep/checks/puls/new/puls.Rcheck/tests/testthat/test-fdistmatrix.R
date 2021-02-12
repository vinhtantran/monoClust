library(fda)

test_that("error for fdistmatrix", {
  expect_error({
    n <- 51
    argvals <- seq(0, 1, len = n)
    #  The true curve values are sine function values with period 1/2
    x <- sin(4 * pi * argvals)
    #  Add independent Gaussian errors with std. dev. 0.2 to the true values
    sigerr <- 0.2
    y <- x + rnorm(x) * sigerr
    #  When we ran this code, we got these values of y (rounded to two
    #  decimals):
    y <- c(0.27,  0.05,  0.58,  0.91,  1.07,  0.98,  0.54,  0.94,  1.13,  0.64,
           0.64,  0.60,  0.24,  0.15, -0.20, -0.63, -0.40, -1.22, -1.11, -0.76,
           -1.11, -0.69, -0.54, -0.50, -0.35, -0.15,  0.27,  0.35,  0.65,  0.75,
           0.75,  0.91,  1.04,  1.04,  1.04,  0.46,  0.30, -0.01, -0.19, -0.42,
           -0.63, -0.78, -1.01, -1.08, -0.91, -0.92, -0.72, -0.84, -0.38, -0.23,
           0.02)
    #  Set up a B-spline basis system of order 4 (piecewise cubic) and with
    #  knots at 0, 0.1, ..., 0.9 and 1.0, and plot the basis functions
    nbasis <- 13
    basisobj <- fda::create.bspline.basis(c(0, 1), nbasis)
    #  Smooth the data, outputting only the functional data object for the
    #  fitted curve.  Note that in this simple case we can supply the basis
    #  object as the "fdParobj" parameter
    ys <- fda::smooth.basis(argvals = argvals, y = y, fdParobj = basisobj)
    fdistmatrix(ys)
  }, "\"fd\" must be of class \"fd\"")
})

test_that("fdistmatrix is correct with usc", {
  expect_equal({
    # Examples taken from fda::Data2fd()
    data(gait)
    # Function only works on two dimensional data
    gait <- gait[, 1:5, 1]
    gaitbasis3 <- create.fourier.basis(nbasis = 5)
    gaitfd3 <- suppressMessages(Data2fd(gait, basisobj = gaitbasis3))

    a <- fdistmatrix(gaitfd3, c(0.2, 0.4), "usc")
    signif(a, 7)
  }, {
    b <- matrix(c(0.000000, 1.508725, 2.702666, 2.959672, 6.555149,
                  1.508725, 0.000000, 2.089777, 3.740105, 7.326595,
                  2.702666, 2.089777, 0.000000, 5.570693, 9.180862,
                  2.959672, 3.740105, 5.570693, 0.000000, 3.643998,
                  6.555149, 7.326595, 9.180862, 3.643998, 0.000000), 5, 5)

    colnames(b) <- rownames(b) <- paste0("boy", 1:5)
    b
  })
})

test_that("fdistmatrix is correct with manual", {
  expect_equal({
    # Examples taken from fda::Data2fd()
    data(gait)
    # Function only works on two dimensional data
    gait <- gait[, 1:5, 1]
    gaitbasis3 <- create.fourier.basis(nbasis = 5)
    gaitfd3 <- suppressMessages(Data2fd(gait, basisobj = gaitbasis3))

    a <- fdistmatrix(gaitfd3, c(0.2, 0.4), "manual")
    signif(a, 7)
  }, {
    a <- matrix(c(0.000000, 1.453237, 2.703014, 2.982000, 6.553343,
                  1.453237, 0.000000, 2.112801, 3.733874, 7.296259,
                  2.703014, 2.112801, 0.000000, 5.603629, 9.191510,
                  2.982000, 3.733874, 5.603629, 0.000000, 3.611041,
                  6.553343, 7.296259, 9.191510, 3.611041, 0.000000), 5, 5)

    colnames(a) <- rownames(a) <- paste0("boy", 1:5)
    a
  })
})

detach("package:fda")
