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

# FCVARoptionUpdates(opt, p, r)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1)

# GetBounds(opt)

opt <- FCVARoptions()
UB_LB_bounds <- GetBounds(opt)

opt <- FCVARoptions()
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
UB_LB_bounds <- GetBounds(opt)


# FCVARestn(x,k,r,opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)



################################################################################
# Specification Functions
################################################################################


# LagSelect(x, kmax, r, order, opt )

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
LagSelect(x, kmax = 3, r = 3, order = 12, opt)

# RankTests(x, k, opt)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
rankTestStats <- RankTests(x, k = 2, opt)


# GetPvalues(q, b, consT, testStat, opt)

# TODO: Make this work.

opt <- FCVARoptions()
pv <- GetPvalues(q = 1, b = 0.4, consT = 0, testStat = 3.84, opt)

opt <- FCVARoptions()
pv <- GetPvalues(q = 1, b = 0.75, consT = 0, testStat = 3.84, opt)



# FCVARbootRank <- function(x, k, opt, r1, r2, B)

# TODO: Test this more extensively.

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
FCVARbootRank_out <- FCVARbootRank(x, k = 2, opt, r1 = 0, r2 = 1, B = 999)



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
MVWNtest_out <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)

set.seed(27)
WN <- rnorm(100)
RW <- cumsum(rnorm(100))
MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
MVWNtest_out <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)


# LMtest(x,q)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
MVWNtest_out <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
LMtest(x = matrix(results$Residuals[, 1]), q = 12)
LMtest(x = results$Residuals[,2, drop = FALSE], q = 12)

set.seed(27)
WN <- rnorm(100)
RW <- cumsum(rnorm(100))
LMtest(x = matrix(WN), q = 10)
LMtest(x = matrix(RW), q = 10)
MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
MVWNtest_out <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)




# Qtest(x,q)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
MVWNtest_out <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
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



# HypoTest(modelUNR, modelR)

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
m1r1 <- FCVARestn(x1, k, r, opt1)
Hdb <- HypoTest(modelUNR = m1, modelR = m1r1)

opt1 <- opt
opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
m1r2 <- FCVARestn(x1, k, r, opt1)
Hbeta1 <- HypoTest(m1, m1r2)

opt1 <- opt
opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
m1r4 <- FCVARestn(x1, k, r, opt1)
Halpha2 <- HypoTest(m1, m1r4)


# FCVARforecast(x, model, NumPeriods)

# TODO: Test this and choose example.


# FCVARboot(x, k, r, optRES, optUNR, B)

# TODO: Test this and choose example.


# GetCharPolyRoots(coeffs, opt, k, r, p)

opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)
cPolyRoots <- GetCharPolyRoots(results$coeffs, opt, k = 2, r = 1, p = 3)


################################################################################
# Auxilliary Functions
################################################################################



################################################################################
# End
################################################################################

