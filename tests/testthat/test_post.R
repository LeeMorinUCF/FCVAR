
context("Postestimation")



test_that("Hypothesis testing results and output are correct", {
  load(file = 'soln_post/HypoTest.RData')
  Hdb_text_soln <- readLines('soln_post/HypoTest_Hdb.txt')
  Hbeta1_text_soln <- readLines('soln_post/HypoTest_Hbeta1.txt')
  Halpha2_text_soln <- readLines('soln_post/HypoTest_Halpha2.txt')

  capture.output(Hdb_test <- HypoTest(modelUNR = m1, modelR = m1r1),
                 file = 'soln_post/temp.txt')
  Hdb_text <- readLines('soln_post/temp.txt')

  expect_equal(Hdb_test, Hdb)
  expect_equal(Hdb_text, Hdb_text_soln)

  capture.output(Hbeta1_test <- HypoTest(m1, m1r2),
                 file = 'soln_post/temp.txt')
  Hbeta1_text <- readLines('soln_post/temp.txt')

  expect_equal(Hbeta1_test, Hbeta1)
  expect_equal(Hbeta1_text, Hbeta1_text_soln)

  capture.output(Halpha2_test <- HypoTest(m1, m1r4),
                 file = 'soln_post/temp.txt')
  Halpha2_text <- readLines('soln_post/temp.txt')

  expect_equal(Halpha2_test, Halpha2)
  expect_equal(Halpha2_text, Halpha2_text_soln)
})


test_that("Bootstrap hypothesis testing results and output are correct", {

  skip('Bootstrap hypothesis test skipped')

  load(file = 'soln_post/FCVARboot_stats.RData')
  FCVARboot_stats_text_soln <- readLines('soln_post/FCVARboot_stats.txt')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  opt$plotRoots <- 0
  optUNR <- opt
  optRES <- opt
  optRES$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
  set.seed(42)
  # FCVARboot_stats <- FCVARboot(x, k = 2, r = 1, optRES, optUNR, B = 999)
  # FCVARboot_stats <- FCVARboot(x[1:50, ], k = 2, r = 1, optRES, optUNR, B = 5)

  capture.output(FCVARboot_stats_test <- FCVARboot(x[1:50, ], k = 2, r = 1, optRES, optUNR, B = 5),
                 file = 'soln_post/temp.txt')
  FCVARboot_stats_text <- readLines('soln_post/temp.txt')

  expect_equal(FCVARboot_stats_test, FCVARboot_stats)
  expect_equal(FCVARboot_stats_text, FCVARboot_stats_text_soln)
})


test_that("forecasts are calculated correctly", {
  load(file = 'soln_post/xf.RData')

  # opt <- FCVARoptions()
  # opt$gridSearch   <- 0 # Disable grid search in optimization.
  # opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  # opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  # opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  # opt1 <- opt
  # opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
  # m1r4 <- FCVARestn(x1, k, r, opt1)
  xf_test <- FCVARforecast(x, m1r4, NumPeriods = 12)

  expect_equal(xf_test, xf)
})
