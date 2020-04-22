
context("Estimation")


test_that("estimation options have correct defaults", {
  load(file = 'soln_estn/opt.RData')
  expect_equal(FCVARoptions(), opt)
})


test_that("estimation options are updated correctly", {
  load(file = 'soln_estn/newOpt.RData')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.

  expect_equal(FCVARoptionUpdates(opt, p = 3, r = 1), newOpt)
})


test_that("default bounds on parameter space are correct", {
  load(file = 'soln_estn/UB_LB_bounds_def.RData')

  opt <- FCVARoptions()

  expect_equal(GetBounds(opt), UB_LB_bounds)
})


test_that("modified bounds on parameter space are correct", {
  load(file = 'soln_estn/UB_LB_bounds_mod.RData')

  opt <- FCVARoptions()
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.

  expect_equal(GetBounds(opt), UB_LB_bounds)
})


test_that("base estimation results are correct", {
  load(file = 'soln_estn/results_m1.RData')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]

  expect_equal(FCVARestn(x, k = 2, r = 1, opt), results)
})


test_that("base estimation printed output is correct", {
  results_text_soln <- readLines('soln_estn/results_m1.txt')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  capture.output(results <- FCVARestn(x, k = 2, r = 1, opt), file = 'soln_estn/temp.txt')
  results_text <- readLines('soln_estn/temp.txt')

  expect_equal(results_text, results_text_soln)
})




