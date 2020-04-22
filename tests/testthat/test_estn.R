
context("Estimation")


test_that("estimation options have correct defaults", {
  load(file = 'soln_estn/opt.RData')
  expect_equal(FCVARoptions(), opt)
})

