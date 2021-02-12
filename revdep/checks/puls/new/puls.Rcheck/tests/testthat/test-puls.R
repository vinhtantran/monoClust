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

test_that("calling puls wrongly creates correct errors", {
  testthat::expect_error({
    PULS()
  },
  "\"toclust.fd\" must be a functional data object (fda::is.fd).",
  fixed = TRUE)

  testthat::expect_error({
    PULS("hello")
  },
  "\"toclust.fd\" must be a functional data object (fda::is.fd).",
  fixed = TRUE)

  testthat::expect_error({
    PULS(yfd$fd, minbucket = 5, minsplit = 2)
  },
  "\"minbucket\" must be less than \"minsplit\".",
  fixed = TRUE)
})

test_that("calling puls correctly in a popular case", {
  skip_on_cran()
  testthat::expect_s3_class({
    PULS4_pam <- PULS(toclust.fd = yfd$fd, intervals = intervals,
                      nclusters = 4, method = "pam")
  },
  "PULS")
})

test_that("calling puls when intervals don't have names", {
  skip_on_cran()
  testthat::expect_output({
    rownames(intervals) <- NULL
    PULS4_pam <- PULS(toclust.fd = yfd$fd, intervals = intervals,
                      nclusters = 4, method = "pam")
    print(PULS4_pam)
  },
  "2) 7 15  885.3640 0.8431711")
})

test_that("calling puls with ward method", {
  skip_on_cran()
  testthat::expect_output({
    PULS4_pam <- PULS(toclust.fd = yfd$fd, intervals = intervals,
                      nclusters = 4, method = "ward")
    print(PULS4_pam)
  },
  "2) Jul 15  885.3640 0.8431711 ")
})
