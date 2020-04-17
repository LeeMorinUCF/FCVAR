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



################################################################################
# Auxilliary Functions
################################################################################



################################################################################
# End
################################################################################

