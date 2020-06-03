
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


# test_that("base estimation results are correct", {
#   load(file = 'soln_estn/results_m1.RData')
#
#   opt <- FCVARoptions()
#   opt$gridSearch   <- 0 # Disable grid search in optimization.
#   opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#   opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#   opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#   x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#
#   expect_equal(FCVARestn(x, k = 2, r = 1, opt), results)
# })
#
#
# test_that("base estimation printed output is correct", {
#   results_text_soln <- readLines('soln_estn/results_m1.txt')
#
#   opt <- FCVARoptions()
#   opt$gridSearch   <- 0 # Disable grid search in optimization.
#   opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#   opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#   opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#   x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#   capture.output(results <- FCVARestn(x, k = 2, r = 1, opt), file = 'soln_estn/temp.txt')
#   results_text <- readLines('soln_estn/temp.txt')
#
#   expect_equal(results_text, results_text_soln)
# })

test_that("unrestricted estimation results and output are correct", {
  load(file = 'soln_estn/results_m1.RData')
  results_text_soln <- readLines('soln_estn/results_m1.txt')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  # opt$db0          <- c(0.6, 0.6) # Set starting values for optimization algorithm.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  # results_test <- FCVARestn(x, k = 2, r = 1, opt)
  capture.output(results_test <- FCVARestn(x, k = 2, r = 1, opt), file = 'soln_estn/temp.txt')
  results_text <- readLines('soln_estn/temp.txt')

  expect_equal(results_test, results)
  expect_equal(results_text, results_text_soln)
})


test_that("restricted estimation results and output are correct", {
  load(file = 'soln_estn/results_m1r124.RData')
  m1r1_text_soln <- readLines('soln_estn/results_m1r1.txt')
  m1r2_text_soln <- readLines('soln_estn/results_m1r2.txt')
  m1r4_text_soln <- readLines('soln_estn/results_m1r4.txt')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]

  opt1 <- opt
  opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
  opt1$r_psi <- 1
  capture.output(m1r1_test <- FCVARestn(x, k = 2, r = 1, opt1),
                 file = 'soln_estn/temp.txt')
  m1r1_text <- readLines('soln_estn/temp.txt')

  expect_equal(m1r1_test, m1r1)
  expect_equal(m1r1_text, m1r1_text_soln)

  opt1 <- opt
  opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
  # opt1$db0 <- c(0.67, 0.67) # Set starting values for optimization algorithm.
  capture.output(m1r2_test <- FCVARestn(x, k = 2, r = 1, opt1),
                 file = 'soln_estn/temp.txt')
  m1r2_text <- readLines('soln_estn/temp.txt')

  expect_equal(m1r2_test, m1r2)
  expect_equal(m1r2_text, m1r2_text_soln)

  opt1 <- opt
  opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
  # opt1$db0 <- c(0.575, 0.575) # Set starting values for optimization algorithm.
  capture.output(m1r4_test <- FCVARestn(x, k = 2, r = 1, opt1),
                 file = 'soln_estn/temp.txt')
  m1r4_text <- readLines('soln_estn/temp.txt')

  expect_equal(m1r4_test, m1r4)
  expect_equal(m1r4_text, m1r4_text_soln)
})

