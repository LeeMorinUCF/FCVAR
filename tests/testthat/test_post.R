
context("Postestimation")


#
# test_that("Lag selection results and output are correct", {
#   load(file = 'soln_estn/results_m1.RData')
#   results_text_soln <- readLines('soln_estn/results_m1.txt')
#
#   opt <- FCVARoptions()
#   opt$gridSearch   <- 0 # Disable grid search in optimization.
#   opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#   opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#   opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#   x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#   # results_test <- FCVARestn(x, k = 2, r = 1, opt)
#   capture.output(results_test <- FCVARestn(x, k = 2, r = 1, opt), file = 'soln_estn/temp.txt')
#   results_text <- readLines('soln_estn/temp.txt')
#
#   expect_equal(results_test, results)
#   expect_equal(results_text, results_text_soln)
# })


