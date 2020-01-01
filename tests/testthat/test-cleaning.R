test_that("Discretizing with percentile cutoff works", {
  x <- seq(from=0, to=100, length.out=101)
  result <- sis_discretize(x, cutoff=0.305)
  expect_equal(sum(result), 70)
})

test_that("Discretizing with absolute cutoff works", {
  x <- seq(from=0, to=100, length.out=101)
  result <- sis_discretize(x, cutoff=17, use_percentile=FALSE)
  expect_equal(sum(result), 83)
})

test_that("Discretizing with NAs works", {
  x <- seq(from=0, to=100, length.out=101)
  x[c(1, 3, 17)] <- NA
  result <- sis_discretize(x, cutoff=0.305)
  expect_identical(result[17], as.numeric(NA))
  expect_equal(result[2], 0)
  expect_equal(result[50],1)
})


test_that("Getting sisters works", {
  phy <- ape::rcoal(15)
  sisters <- sis_get_sisters(phy, ncores=1)
  expect_equal(class(sisters), "data.frame")
})
