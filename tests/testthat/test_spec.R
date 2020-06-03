
context("Specification")



test_that("Lag selection results and output are correct", {
  load(file = 'soln_spec/LagSelectStats.RData')
  LagSelectStats_text_soln <- readLines('soln_spec/LagSelectStats.txt')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  # LagSelectStats <- LagSelect(x, kmax = 3, r = 3, order = 12, opt)

  capture.output(LagSelectStats_test <- LagSelect(x, kmax = 3, r = 3, order = 12, opt),
                 file = 'soln_spec/temp.txt')
  LagSelectStats_text <- readLines('soln_spec/temp.txt')

  expect_equal(LagSelectStats_test, LagSelectStats)
  expect_equal(LagSelectStats_text, LagSelectStats_text_soln)
})


test_that("Rank Testing results and output are correct", {
  load(file = 'soln_spec/rankTestStats.RData')
  rankTestStats_text_soln <- readLines('soln_spec/rankTestStats.txt')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  # rankTestStats <- RankTests(x, k = 2, opt)

  capture.output(rankTestStats_test <- RankTests(x, k = 2, opt),
                 file = 'soln_spec/temp.txt')
  rankTestStats_text <- readLines('soln_spec/temp.txt')

  expect_equal(rankTestStats_test, rankTestStats)
  expect_equal(rankTestStats_text, rankTestStats_text_soln)
})



test_that("Bootstrap Rank Testing results and output are correct", {
  load(file = 'soln_spec/FCVARbootRank_stats.RData')
  FCVARbootRank_stats_text_soln <- readLines('soln_spec/FCVARbootRank_stats.txt')

  opt <- FCVARoptions()
  opt$gridSearch   <- 0 # Disable grid search in optimization.
  opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
  opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
  opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  opt$plotRoots <- 0
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  set.seed(42)
  # FCVARbootRank_stats <- FCVARbootRank(x, k = 2, opt, r1 = 0, r2 = 1, B = 999)
  # FCVARbootRank_stats <- FCVARbootRank(x[1:50, ], k = 2, opt, r1 = 0, r2 = 1, B = 5)

  capture.output(FCVARbootRank_stats_test <- FCVARbootRank(x[1:50, ], k = 2, opt, r1 = 0, r2 = 1, B = 5),
                 file = 'soln_spec/temp.txt')
  FCVARbootRank_stats_text <- readLines('soln_spec/temp.txt')

  expect_equal(FCVARbootRank_stats_test, FCVARbootRank_stats)
  expect_equal(FCVARbootRank_stats_text, FCVARbootRank_stats_text_soln)
})


