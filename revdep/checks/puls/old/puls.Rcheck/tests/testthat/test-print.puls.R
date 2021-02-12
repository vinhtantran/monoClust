library(fda)

# Build a simple fd object from already smoothed smoothed_arctic
data(smoothed_arctic)
NBASIS <- 300
NORDER <- 4
y <- t(as.matrix(smoothed_arctic[, -1]))
splinebasis <- create.bspline.basis(rangeval = c(1, 365),
                                    nbasis = NBASIS,
                                    norder = NORDER)
fdParobj <- fdPar(fdobj = splinebasis,
                  Lfdobj = 2,
                  # No need for any more smoothing
                  lambda = .000001)
yfd <- smooth.basis(argvals = 1:365, y = y, fdParobj = fdParobj)

Jan <- c(1, 31); Feb <- c(31, 59); Mar <- c(59, 90)
Apr <- c(90, 120); May <- c(120, 151); Jun <- c(151, 181)
Jul <- c(181, 212); Aug <- c(212, 243); Sep <- c(243, 273)
Oct <- c(273, 304); Nov <- c(304, 334); Dec <- c(334, 365)

intervals <-
  rbind(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec)

PULS4_pam <- PULS(toclust.fd = yfd$fd, intervals = intervals,
                  nclusters = 4, method = "pam")


test_that("print.puls creates correct output", {
  skip_on_cran()

  testthat::expect_output(
    print(PULS4_pam),
    paste("Note: One or more of the splits chosen had an alternative split",
          "that reduced inertia by the same amount. See \"alt\" column of",
          "\"frame\" object for details.")
  )

  testthat::expect_output(
    print(PULS4_pam),
    "4) Aug 8  311.7792   *"
  )
})

detach("package:fda")
