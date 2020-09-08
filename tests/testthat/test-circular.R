test_that("circular data works", {
  expect_equal(MonoClust(data.frame(wind.sensit.bin.2007), cir.var = 3, nclusters = 2), 2)
})
