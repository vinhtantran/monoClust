test_that("correct mean fd calculation", {
  expect_equal({
    bspl4 <- fda::create.bspline.basis(nbasis = 4)
    parab4.5 <- fda::fd(c(3, -1, -1, 3) / 3, bspl4)
    mean_fd(parab4.5)
  }, 1)
})

test_that("incorrect data are transferred", {
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
    mean_fd(ys)
  }, "'x' is not of class 'fd'")
})
