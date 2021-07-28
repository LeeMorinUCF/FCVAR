
context("Estimation")


test_that("estimation options have correct defaults", {
  load(file = 'soln_estn/opt.RData')
  expect_equal(FCVARoptions(), opt)
})


test_that("estimation options are updated correctly", {
  load(file = 'soln_estn/newOpt.RData')

  opt <- FCVARoptions(
    gridSearch   = 0, # Disable grid search in optimization.
    dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
    dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
    constrained  = 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  )

  expect_equal(FCVARoptionUpdates(opt, p = 3, r = 1), newOpt)
})


test_that("default bounds on parameter space are correct", {
  load(file = 'soln_estn/UB_LB_bounds_def.RData')

  opt <- FCVARoptions()

  expect_equal(GetBounds(opt), UB_LB_bounds)
})


test_that("modified bounds on parameter space are correct", {
  load(file = 'soln_estn/UB_LB_bounds_mod.RData')

  opt <- FCVARoptions(
    dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
    dbMax        = c(2.00, 2.00) # Set upper bound for d,b.
  )

  expect_equal(GetBounds(opt), UB_LB_bounds)
})


test_that("unrestricted estimation results and output are correct", {
  load(file = 'soln_estn/results_m1.RData')
  results_text_soln <- readLines('soln_estn/results_m1.txt')

  opt <- FCVARoptions(
    unc_optim_control = list(maxit = 1000,    # Reduce tolerance.
                             reltol = 1e-10),
    con_optim_control = list(maxit = 1000,    # Reduce tolerance.
                             pgtol = 1e-10),
    hess_delta = 10^(-3),  # Stable for testing on 32-bit, 386 architecture.
    gridSearch   = 0, # Disable grid search in optimization.
    dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
    dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
    constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
    plotRoots    = 0 # Don't create plots for tests.
  )
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
  # results_test <- FCVARestn(x, k = 2, r = 1, opt)
  capture.output(results_test <- FCVARestn(x, k = 2, r = 1, opt), file = 'soln_estn/temp.txt')
  results_text <- readLines('soln_estn/temp.txt')

  # Note that the inverse Hessian matrix has cross-platform differences in
  # mean relative difference of 1.968147e-08, which is not a problem.
  # results$NegInvHessian <- round(results$NegInvHessian, 6)
  # results_test$NegInvHessian <- round(results_test$NegInvHessian, 6)
  # # Similarly for other estimates.
  # results$SE$db <- round(results$SE$db, 6)
  # results_test$SE$db <- round(results_test$SE$db, 6)
  # results$SE$muHat <- round(results$SE$muHat, 6)
  # results_test$SE$muHat <- round(results_test$SE$muHat, 6)
  # results$SE$alphaHat <- round(results$SE$alphaHat, 6)
  # results_test$SE$alphaHat <- round(results_test$SE$alphaHat, 6)
  # # results$SE$gammaHat <- round(results$SE$gammaHat, 6)
  # # results_test$SE$gammaHat <- round(results_test$SE$gammaHat, 6)
  # Conduct separate tests on Hessian and standard errors.
  results_NegInvHessian <- results$NegInvHessian
  results$NegInvHessian <- NULL
  results_test_NegInvHessian <- results_test$NegInvHessian
  results_test$NegInvHessian <- NULL
  results_SE <- results$SE
  results$SE <- NULL
  results_test_SE <- results_test$SE
  results_test$SE <- NULL
  # Accuracy can only be so high on 32-bit platforms.

  expect_equal(results_test, results)
  expect_equal(results_test_NegInvHessian, results_NegInvHessian, tolerance = 10^(-5))
  expect_equal(results_test_SE, results_SE, tolerance = 10^(-6))
  expect_equal(results_text, results_text_soln)
})


test_that("restricted estimation results and output are correct", {
  load(file = 'soln_estn/results_m1r124.RData')
  m1r1_text_soln <- readLines('soln_estn/results_m1r1.txt')
  m1r2_text_soln <- readLines('soln_estn/results_m1r2.txt')
  m1r4_text_soln <- readLines('soln_estn/results_m1r4.txt')

  opt <- FCVARoptions(
    unc_optim_control = list(maxit = 1000,    # Reduce tolerance.
                             reltol = 1e-10),
    con_optim_control = list(maxit = 1000,    # Reduce tolerance.
                             pgtol = 1e-10),
    hess_delta = 10^(-3),  # Stable for testing on 32-bit, 386 architecture.
    gridSearch   = 0, # Disable grid search in optimization.
    dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
    dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
    constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
    plotRoots    = 0 # Don't create plots for tests.
  )
  x <- votingJNP2014[, c("lib", "ir_can", "un_can")]

  opt1 <- opt
  opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
  opt1$r_psi <- 1
  capture.output(m1r1_test <- FCVARestn(x, k = 2, r = 1, opt1),
                 file = 'soln_estn/temp.txt')
  m1r1_text <- readLines('soln_estn/temp.txt')

  # Note that the inverse Hessian matrix has cross-platform differences in
  # mean relative difference on the order of 1.0e-08, which is not a problem.
  # m1r1$NegInvHessian <- round(m1r1$NegInvHessian, 6)
  # m1r1_test$NegInvHessian <- round(m1r1_test$NegInvHessian, 6)
  # # Similarly for other estimates.
  # m1r1$SE$db <- round(m1r1$SE$db, 6)
  # m1r1_test$SE$db <- round(m1r1_test$SE$db, 6)
  # m1r1$SE$muHat <- round(m1r1$SE$muHat, 6)
  # m1r1_test$SE$muHat <- round(m1r1_test$SE$muHat, 6)
  # m1r1$SE$alphaHat <- round(m1r1$SE$alphaHat, 6)
  # m1r1_test$SE$alphaHat <- round(m1r1_test$SE$alphaHat, 6)
  # # m1r1$SE$gammaHat <- round(m1r1$SE$gammaHat, 6)
  # # m1r1_test$SE$gammaHat <- round(m1r1_test$SE$gammaHat, 6)
  # Conduct separate tests on Hessian and standard errors.
  m1r1_NegInvHessian <- m1r1$NegInvHessian
  m1r1$NegInvHessian <- NULL
  m1r1_test_NegInvHessian <- m1r1_test$NegInvHessian
  m1r1_test$NegInvHessian <- NULL
  m1r1_SE <- m1r1$SE
  m1r1$SE <- NULL
  m1r1_test_SE <- m1r1_test$SE
  m1r1_test$SE <- NULL
  # Accuracy can only be so high on 32-bit platforms.

  expect_equal(m1r1_test, m1r1)
  expect_equal(m1r1_test_NegInvHessian, m1r1_NegInvHessian, tolerance = 10^(-5))
  expect_equal(m1r1_test_SE, m1r1_SE, tolerance = 10^(-6))
  expect_equal(m1r1_text, m1r1_text_soln)

  opt1 <- opt
  opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
  # opt1$db0 <- c(0.67, 0.67) # Set starting values for optimization algorithm.
  capture.output(m1r2_test <- FCVARestn(x, k = 2, r = 1, opt1),
                 file = 'soln_estn/temp.txt')
  m1r2_text <- readLines('soln_estn/temp.txt')

  # Note that the inverse Hessian matrix has cross-platform differences in
  # mean relative difference on the order of 1.0e-08, which is not a problem.
  # m1r2$NegInvHessian <- round(m1r2$NegInvHessian, 6)
  # m1r2_test$NegInvHessian <- round(m1r2_test$NegInvHessian, 6)
  # # Similarly for other estimates.
  # m1r2$SE$db <- round(m1r2$SE$db, 6)
  # m1r2_test$SE$db <- round(m1r2_test$SE$db, 6)
  # m1r2$SE$muHat <- round(m1r2$SE$muHat, 6)
  # m1r2_test$SE$muHat <- round(m1r2_test$SE$muHat, 6)
  # m1r2$SE$alphaHat <- round(m1r2$SE$alphaHat, 6)
  # m1r2_test$SE$alphaHat <- round(m1r2_test$SE$alphaHat, 6)
  # # m1r2$SE$gammaHat <- round(m1r2$SE$gammaHat, 6)
  # # m1r2_test$SE$gammaHat <- round(m1r2_test$SE$gammaHat, 6)
  # Conduct separate tests on Hessian and standard errors.
  m1r2_NegInvHessian <- m1r2$NegInvHessian
  m1r2$NegInvHessian <- NULL
  m1r2_test_NegInvHessian <- m1r2_test$NegInvHessian
  m1r2_test$NegInvHessian <- NULL
  m1r2_SE <- m1r2$SE
  m1r2$SE <- NULL
  m1r2_test_SE <- m1r2_test$SE
  m1r2_test$SE <- NULL
  # Accuracy can only be so high on 32-bit platforms.

  expect_equal(m1r2_test, m1r2)
  expect_equal(m1r2_test_NegInvHessian, m1r2_NegInvHessian, tolerance = 10^(-5))
  expect_equal(m1r2_test_SE, m1r2_SE, tolerance = 10^(-6))
  expect_equal(m1r2_text, m1r2_text_soln)

  # expect_equal(m1r2_test, m1r2)
  #
  # # Similarly, some numbers in SE different in last decimal place.
  # expect_equal(m1r2_text, m1r2_text_soln)
  # # expect_equal(m1r2_text[1:85], m1r2_text_soln[1:85])
  # # expect_equal(m1r2_text[87:103], m1r2_text_soln[87:103])
  # # expect_equal(m1r2_text[105:113], m1r2_text_soln[105:113])

  opt1 <- opt
  opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
  # opt1$db0 <- c(0.575, 0.575) # Set starting values for optimization algorithm.
  capture.output(m1r4_test <- FCVARestn(x, k = 2, r = 1, opt1),
                 file = 'soln_estn/temp.txt')
  m1r4_text <- readLines('soln_estn/temp.txt')

  # Note that the inverse Hessian matrix has cross-platform differences in
  # mean relative difference on the order of 1.0e-08, which is not a problem.
  # m1r4$NegInvHessian <- round(m1r4$NegInvHessian, 6)
  # m1r4_test$NegInvHessian <- round(m1r4_test$NegInvHessian, 6)
  # # Similarly for other estimates.
  # m1r4$SE$db <- round(m1r4$SE$db, 6)
  # m1r4_test$SE$db <- round(m1r4_test$SE$db, 6)
  # m1r4$SE$muHat <- round(m1r4$SE$muHat, 6)
  # m1r4_test$SE$muHat <- round(m1r4_test$SE$muHat, 6)
  # m1r4$SE$alphaHat <- round(m1r4$SE$alphaHat, 6)
  # m1r4_test$SE$alphaHat <- round(m1r4_test$SE$alphaHat, 6)
  # # m1r4$SE$gammaHat <- round(m1r4$SE$gammaHat, 6)
  # # m1r4_test$SE$gammaHat <- round(m1r4_test$SE$gammaHat, 6)
  # Conduct separate tests on Hessian and standard errors.
  m1r4_NegInvHessian <- m1r4$NegInvHessian
  m1r4$NegInvHessian <- NULL
  m1r4_test_NegInvHessian <- m1r4_test$NegInvHessian
  m1r4_test$NegInvHessian <- NULL
  m1r4_SE <- m1r4$SE
  m1r4$SE <- NULL
  m1r4_test_SE <- m1r4_test$SE
  m1r4_test$SE <- NULL
  # Accuracy can only be so high on 32-bit platforms.

  expect_equal(m1r4_test, m1r4)
  expect_equal(m1r4_test_NegInvHessian, m1r4_NegInvHessian, tolerance = 10^(-5))
  expect_equal(m1r4_test_SE, m1r4_SE, tolerance = 10^(-6))
  expect_equal(m1r4_text, m1r4_text_soln)

  # expect_equal(m1r4_test, m1r4)
  #
  # # Similarly, some numbers in SE different in last decimal place.
  # expect_equal(m1r4_text, m1r4_text_soln)
  # # expect_equal(m1r4_text[1:85], m1r4_text_soln[1:85])
  # # expect_equal(m1r4_text[87:103], m1r4_text_soln[87:103])
  # # expect_equal(m1r4_text[105:113], m1r4_text_soln[105:113])
})

