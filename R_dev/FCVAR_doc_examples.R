################################################################################
#
# Examples of Functions in the FCVAR Model
#
#
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
#
# April 17, 2020
#
################################################################################
#
# This is a series of examples of the usage of each function in the FCVAR package.
# It is the list of examples used in the documentation for the functions.
#
################################################################################



################################################################################
# Estimation Functions
################################################################################

# FCVARoptions

opt <- FCVARoptions()
# save(opt, list = c('opt'), file = 'tests/testthat/soln_estn/opt.RData')


# FCVARoptionUpdates(opt, p, r)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1)
# save(newOpt, list = c('newOpt'), file = 'tests/testthat/soln_estn/newOpt.RData')


# GetBounds(opt)

opt <- FCVARoptions()
UB_LB_bounds <- GetBounds(opt)
# save(UB_LB_bounds, list = c('UB_LB_bounds'), file = 'tests/testthat/soln_estn/UB_LB_bounds_def.RData')

opt <- FCVARoptions()
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
UB_LB_bounds <- GetBounds(opt)
# save(UB_LB_bounds, list = c('UB_LB_bounds'), file = 'tests/testthat/soln_estn/UB_LB_bounds_mod.RData')


# FCVARestn(x,k,r,opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
# save(results, list = c('results'), file = 'tests/testthat/soln_estn/results_m1.RData')
# capture.output(results <- FCVARestn(x, k = 2, r = 1, opt), file = 'tests/testthat/soln_estn/results_m1.txt')

opt1 <- opt
opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt1$r_psi <- 1
m1r1 <- FCVARestn(x, k = 2, r = 1, opt1)
# capture.output(m1r1 <- FCVARestn(x, k = 2, r = 1, opt1),
#                file = 'tests/testthat/soln_estn/results_m1r1.txt')

opt1 <- opt
opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
m1r2 <- FCVARestn(x, k = 2, r = 1, opt1)
# capture.output(m1r2 <- FCVARestn(x, k = 2, r = 1, opt1),
#                file = 'tests/testthat/soln_estn/results_m1r2.txt')

opt1 <- opt
opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
m1r4 <- FCVARestn(x, k = 2, r = 1, opt1)
# capture.output(m1r4 <- FCVARestn(x, k = 2, r = 1, opt1),
#                file = 'tests/testthat/soln_estn/results_m1r4.txt')

# save(m1r1, m1r2, m1r4,
#      list = c('m1r1', 'm1r2', 'm1r4'),
#      file = 'tests/testthat/soln_estn/results_m1r124.RData')




################################################################################
# Specification Functions
################################################################################


# FCVARlagSelect(x, kmax, r, order, opt )

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
FCVARlagSelectStats <- FCVARlagSelect(x, kmax = 3, r = 3, order = 12, opt)
# capture.output(FCVARlagSelectStats <- FCVARlagSelect(x, kmax = 3, r = 3, order = 12, opt), file = 'tests/testthat/soln_spec/FCVARlagSelectStats.txt')
# save(FCVARlagSelectStats, list = c('FCVARlagSelectStats'), file = 'tests/testthat/soln_spec/FCVARlagSelectStats.RData')


# print.FCVARlagSelect(stats, kmax, r, p, order, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
FCVARlagSelectStats <- FCVARlagSelect(x, kmax = 3, r = 3, order = 12, opt)
print.FCVARlagSelect(stats = FCVARlagSelectStats, kmax = 3, r = 3, p = 3, order, opt)


# FCVARrankTests(x, k, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
rankTestStats <- FCVARrankTests(x, k = 2, opt)
# capture.output(rankTestStats <- FCVARrankTests(x, k = 2, opt), file = 'tests/testthat/soln_spec/rankTestStats.txt')
# save(rankTestStats, list = c('rankTestStats'), file = 'tests/testthat/soln_spec/rankTestStats.RData')


# print.FCVARrankTests <- function(stats, k, p, T, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
rankTestStats <- FCVARrankTests(x, k = 2, opt)
print.FCVARrankTests(stats = rankTestStats, k = 2, p = ncol(x), T = nrow(x), opt)


# FCVARbootRank <- function(x, k, opt, r1, r2, B)

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
# capture.output(FCVARbootRank_stats <- FCVARbootRank(x[1:50, ], k = 2, opt, r1 = 0, r2 = 1, B = 5), file = 'tests/testthat/soln_spec/FCVARbootRank_stats.txt')
# save(FCVARbootRank_stats, list = c('FCVARbootRank_stats'), file = 'tests/testthat/soln_spec/FCVARbootRank_stats.RData')



################################################################################
# Postestimation Functions
################################################################################


# MVWNtest(x, maxlag, printResults)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)

set.seed(27)
WN <- rnorm(100)
RW <- cumsum(rnorm(100))
MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
MVWNtest_stats <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)


# print.MVWNtest(stats, maxlag)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
print.MVWNtest(stats = MVWNtest_stats, maxlag = 12)

set.seed(27)
WN <- rnorm(100)
RW <- cumsum(rnorm(100))
MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
MVWNtest_stats <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)
print.MVWNtest(stats = MVWNtest_stats, maxlag = 12)




# LMtest(x,q)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
LMtest(x = matrix(results$Residuals[, 1]), q = 12)
LMtest(x = results$Residuals[,2, drop = FALSE], q = 12)

set.seed(27)
WN <- rnorm(100)
RW <- cumsum(rnorm(100))
LMtest(x = matrix(WN), q = 10)
LMtest(x = matrix(RW), q = 10)
MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
MVWNtest_stats <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)




# Qtest(x,q)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
Qtest(x = results$Residuals, maxlag = 12)
Qtest(x = matrix(results$Residuals[, 1]), maxlag = 12)
Qtest(x = results$Residuals[,2, drop = FALSE], maxlag = 12)

set.seed(27)
WN <- rnorm(100)
RW <- cumsum(rnorm(100))
MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
Qtest(x = MVWN_x, maxlag = 10)
Qtest(x = matrix(WN), maxlag = 10)
Qtest(x = matrix(RW), maxlag = 10)



# FCVARhypoTest(modelUNR, modelR)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
m1 <- FCVARestn(x, k = 2, r = 1, opt)
opt1 <- opt
opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt1$r_psi <- 1
m1r1 <- FCVARestn(x, k = 2, r = 1, opt1)
Hdb <- FCVARhypoTest(modelUNR = m1, modelR = m1r1)

opt1 <- opt
opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
m1r2 <- FCVARestn(x, k = 2, r = 1, opt1)
Hbeta1 <- FCVARhypoTest(m1, m1r2)

opt1 <- opt
opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
m1r4 <- FCVARestn(x, k = 2, r = 1, opt1)
Halpha2 <- FCVARhypoTest(m1, m1r4)

# capture.output(Hdb <- FCVARhypoTest(modelUNR = m1, modelR = m1r1),
#                file = 'tests/testthat/soln_post/FCVARhypoTest_Hdb.txt')
# capture.output(Hbeta1 <- FCVARhypoTest(m1, m1r2),
#                file = 'tests/testthat/soln_post/FCVARhypoTest_Hbeta1.txt')
# capture.output(Halpha2 <- FCVARhypoTest(m1, m1r4),
#                file = 'tests/testthat/soln_post/FCVARhypoTest_Halpha2.txt')
# save(m1, m1r1, m1r2, m1r4, Hdb, Hbeta1, Halpha2,
#      list = c('m1', 'm1r1', 'm1r2', 'm1r4', 'Hdb', 'Hbeta1', 'Halpha2'),
#      file = 'tests/testthat/soln_post/FCVARhypoTest.RData')


# Test warnings and errors from corner cases of hypothesis tests.
# Without errors:
LR_test <- FCVARhypoTest(modelUNR = list(like = 3, fp = 22),
                         modelR = list(like = 0, fp = 20))
# Same likelihood value:
LR_test <- FCVARhypoTest(modelUNR = list(like = 3, fp = 22),
                         modelR = list(like = 3, fp = 20))
# Unordered likelihood value:
LR_test <- FCVARhypoTest(modelUNR = list(like = 3, fp = 22),
                         modelR = list(like = 6, fp = 20))
# Unordered number of free parameters:
LR_test <- FCVARhypoTest(modelUNR = list(like = 3, fp = 18),
                         modelR = list(like = 0, fp = 20))


# FCVARboot(x, k, r, optRES, optUNR, B)

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
FCVARboot_stats <- FCVARboot(x, k = 2, r = 1, optRES, optUNR, B = 999)
# FCVARboot_stats <- FCVARboot(x[1:50, ], k = 2, r = 1, optRES, optUNR, B = 5)
# capture.output(FCVARboot_stats <- FCVARboot(x[1:50, ], k = 2, r = 1, optRES, optUNR, B = 5),
#                file = 'tests/testthat/soln_post/FCVARboot_stats.txt')
# save(FCVARboot_stats, list = c('FCVARboot_stats'),
#      file = 'tests/testthat/soln_post/FCVARboot_stats.RData')


# FCVARforecast(x, model, NumPeriods)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
opt1 <- opt
opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
m1r4 <- FCVARestn(x, k = 2, r = 1, opt1)
xf <- FCVARforecast(x, m1r4, NumPeriods = 12)
# save(xf, m1r4, list = c('xf', 'm1r4'), file = 'tests/testthat/soln_post/xf.RData')


# GetCharPolyRoots(coeffs, opt, k, r, p)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
cPolyRoots <- GetCharPolyRoots(results$coeffs, opt, k = 2, r = 1, p = 3)


# print.GetCharPolyRoots(cPolyRoots)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
cPolyRoots <- GetCharPolyRoots(results$coeffs, opt, k = 2, r = 1, p = 3)
print.GetCharPolyRoots(cPolyRoots)
plot.GetCharPolyRoots(cPolyRoots, b = results$coeffs$db[2])

# plot.GetCharPolyRoots(cPolyRoots)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
cPolyRoots <- GetCharPolyRoots(results$coeffs, opt, k = 2, r = 1, p = 3)
print.GetCharPolyRoots(cPolyRoots)
plot.GetCharPolyRoots(cPolyRoots, b= results$coeffs$db[2])


################################################################################
# Auxilliary Functions
################################################################################

# FCVARsim(x, model, NumPeriods)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
x_sim <- FCVARsim(x, results, T_sim = 100)


# FCVARsimBS(data, model, NumPeriods)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
xBS <- FCVARsimBS(x[1:10, ], results, NumPeriods = 100)


# GetParams <- function(x, k, r, db, opt)

# TODO: Revise this example after the global variable is removed.

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)

# estimatesTEMP <- get('estimatesTEMP', envir = .GlobalEnv)
#
# y <- x - matrix(1, nrow = T+opt$N, ncol = p) %*% diag(results$coeffs$muHat)
#
# estimates <- GetParams(y, k = 2, r = 1, db = results$coeffs$db, opt)


opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
opt$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt$r_psi <- 1
results <- FCVARestn(x, k = 2, r = 1, opt)

y <- x - matrix(1, nrow = T+opt$N, ncol = p) %*% diag(results$coeffs$muHat)

estimates <- GetParams(y, k = 2, r = 1, db = results$coeffs$db, opt)



# FCVARlikeGrid(x, k, r, opt)

# Restrict equality of fractional parameters.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.1 # Coarser grid for plotting example.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 1 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$progress     <- 1 # Show progress bar update on each value of b.
# opt$progress     <- 2 # Show progress report on each value of b.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
# plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt, main = 'default')
plot(x = likeGrid_params)
# Add a point for the optimum.
points(likeGrid_params$bGrid_orig[which.max(likeGrid_params$like)],
       max(likeGrid_params$like),
       col = 'red', pch = 16)
# Compare likelihood values.
# cbind(likeGrid_params$dGrid, likeGrid_params$bGrid, likeGrid_params$like)
likeGrid_params$dbHatStar
# cbind(likeGrid_params$dGrid_orig, likeGrid_params$bGrid_orig, likeGrid_params$like)
likeGrid_params$local_max
likeGrid_params$max_like


# Linear restriction on fractional parameters.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.1 # Coarser grid for plotting example.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
# Impose linear restriction on d and b:
opt$R_psi        <- matrix(c(2, -1), nrow = 1, ncol = 2)
opt$r_psi        <- 0.5
opt$progress     <- 1 # Show progress bar update on each value of b.
# opt$progress     <- 2 # Show progress report on each value of b.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
# plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt, main = 'default')
plot(x = likeGrid_params)
# Add a point for the optimum.
points(likeGrid_params$bGrid[which.max(likeGrid_params$like)],
       max(likeGrid_params$like),
       col = 'red', pch = 16)
points(likeGrid_params$bGrid[which(likeGrid_params$like == likeGrid_params$max_like)],
       max(likeGrid_params$max_like),
       col = 'red', pch = 16, cex = 2)
# Compare likelihood values.
# cbind(likeGrid_params$dGrid, likeGrid_params$bGrid, likeGrid_params$like)
likeGrid_params$dbHatStar
likeGrid_params$max_like
# cbind(likeGrid_params$dGrid_orig, likeGrid_params$bGrid_orig, likeGrid_params$like)
likeGrid_params$local_max



# Constrained 2-dimensional optimization.
# Impose restriction dbMax >= d >= b >= dbMin.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.1 # Coarser grid for plotting example.
opt$dbStep2D     <- 0.2 # Coarser grid for plotting example.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 1 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$progress     <- 1 # Show progress bar update on each value of b.
# opt$progress     <- 2 # Show progress report on each value of b.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
# max(likeGrid_params$like, na.rm = TRUE)
likeGrid_params$dbHatStar
likeGrid_params$max_like

# Unconstrained 2-dimensional optimization.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.01 # Coarser grid for plotting example.
opt$dbStep2D     <- 0.2 # Coarser grid for plotting example.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$progress     <- 1 # Show progress bar update on each value of b.
# opt$progress     <- 2 # Show progress report on each value of b.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)




# plot.FCVAR_grid(likeGrid_params, ...)
# Previously:
# plot.FCVARlikeGrid(likeGrid_params, k, r, opt, file, file_ext, main)

# Restrict equality of fractional parameters.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.1 # Coarser grid for plotting example.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 1 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$progress     <- 2 # Show progress report on each value of b.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
# plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt, main = 'default')
plot(x = likeGrid_params)

# FCVARlikeMu(mu, y, db, k, r, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
FCVARlikeMu(mu = results$coeffs$muHat, y = x, db = results$coeffs$db, k = 2, r = 1, opt)


# FCVARlike(params, x, k, r, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
FCVARlike(c(results$coeffs$db, results$coeffs$muHat), x, k = 2, r = 1, opt)


# FCVARlikeFull(x, k, r, coeffs, beta, rho, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
FCVARlikeFull(x, k = 2, r = 1, coeffs = results$coeffs,
              beta = results$coeffs$betaHat, rho = results$coeffs$rhoHat, opt)


# TransformData(x, k, db, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
Z_array <- TransformData(x, k = 2, db = results$coeffs$db, opt)


# GetResiduals(x, k, r, coeffs, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
epsilon <- GetResiduals(x, k = 2, r = 1, coeffs = results$coeffs, opt)


# Lbk(x, b, k)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
Lbkx <- Lbk(x, b = results$coeffs$db[2], k = 2)


# FracDiff(x, d)

set.seed(42)
WN <- matrix(rnorm(200), nrow = 100, ncol = 2)
MVWNtest_stats <- MVWNtest(x = WN, maxlag = 10, printResults = 1)
x <- FracDiff(x = WN, d = - 0.5)
MVWNtest_stats <- MVWNtest(x = x, maxlag = 10, printResults = 1)
WN_x_d <- FracDiff(x, d = 0.5)
MVWNtest_stats <- MVWNtest(x = WN_x_d, maxlag = 10, printResults = 1)

# Compare with results from other packages.
x <- matrix(cumsum(rnorm(100)), nrow = 100, ncol = 1)
d <- 0.73
summary(fracdiff::diffseries(x, d) - FCVAR::FracDiff(x - mean(x), d))
summary(LongMemoryTS::fdiff(x, d) - FCVAR::FracDiff(x, d))


# GetFreeParams(k, r, p, opt, rankJ)

opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = NULL)

opt <- FCVARoptions()
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = 4)

opt <- FCVARoptions()
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = 4)


# FCVARhess(x, k, r, coeffs, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
hessian <- FCVARhess(x, k = 2, r = 1, coeffs = results$coeffs, opt)


# SEmat2vecU

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
params <- SEmat2vecU(coeffs = results$coeffs, k = 2, r = 1, p = 3, opt)
coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )

params <- matrix(seq(25))
coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
params <- SEmat2vecU(coeffs = coeffs, k = 2, r = 1, p = 3, opt)


# SEvec2matU(param, k, r, p, opt )

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
params <- SEmat2vecU(coeffs = results$coeffs, k = 2, r = 1, p = 3, opt)
coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )

params <- matrix(seq(25))
coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
params <- SEmat2vecU(coeffs = coeffs, k = 2, r = 1, p = 3, opt)


# GetRestrictedParams(beta0, S00, S01, S11, T, p, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
opt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
Z_array <- TransformData(x, k = 2, db = results$coeffs$db, opt)
Z0 <- Z_array$Z0
Z1 <- Z_array$Z1
Z2 <- Z_array$Z2
Z3 <- Z_array$Z3

T <- nrow(Z0)
p <- ncol(Z0)
p1 <- ncol(Z1)

Z0hat <- Z0
Z1hat <- Z1
Z2hat <- Z2

R0 <- Z0hat - Z2hat %*% ( solve(t(Z2hat) %*% Z2hat) %*% t(Z2hat) %*% Z0hat )
R1 <- Z1hat - Z2hat %*% ( solve(t(Z2hat) %*% Z2hat) %*% t(Z2hat) %*% Z1hat )

S00 <- t(R0) %*% R0/T
S01 <- t(R0) %*% R1/T
S10 <- t(R1) %*% R0/T
S11 <- t(R1) %*% R1/T

eig_out <- eigen( solve(S11) %*% S10 %*% solve(S00) %*% S01 )
D <- eig_out$values # Only the vector, not the diagonal matrix.
V1 <- eig_out$vectors
V2 <- t( cbind(t(V1), D)[order(D), ] )
V <- V2[1:p1, , drop = FALSE]
betaStar <- V[ 1:p1, seq(p1, p1-r+1, by = -1), drop = FALSE]

switched_mats <- GetRestrictedParams(betaStar, S00, S01, S11, T, p, opt)


# Simpler version: Just give the numbers:
betaStar <- matrix(c(-0.95335616, -0.07345676, -0.29277318), nrow = 3)
S00 <- matrix(c(0.0302086527,  0.001308664,  0.0008200923,
                0.0013086640,  0.821417610, -0.1104617893,
                0.0008200923, -0.110461789,  0.0272861128), nrow = 3)
S01 <- matrix(c(-0.0047314320, -0.04488533,  0.006336798,
                 0.0026708007,  0.17463884, -0.069006455,
                -0.0003414163, -0.07110324,  0.022830494), nrow = 3, byrow = TRUE)
S11 <- matrix(c( 0.061355941,  -0.4109969,  -0.007468716,
                -0.410996895,  70.6110313, -15.865097810,
                -0.007468716, -15.8650978,   3.992435799), nrow = 3)
switched_mats <- GetRestrictedParams(betaStar, S00, S01, S11, T = 316, p = 3, opt)



# Not worth it:
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# opt$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# Z_array <- TransformData(x, k = 2, db = results$coeffs$db, opt)



################################################################################
# End
################################################################################

