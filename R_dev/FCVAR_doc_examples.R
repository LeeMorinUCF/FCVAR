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

# FCVARsimBS(data, model, NumPeriods)

# TODO: Test this and choose example.


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
results <- FCVARestn(x1, k, r, opt)

y <- x - matrix(1, nrow = T+opt$N, ncol = p) %*% diag(results$coeffs$muHat)

estimates <- GetParams(y, k = 2, r = 1, db = results$coeffs$db, opt)



# LikeGridSearch(x, k, r, opt)

# TODO: Test this and choose example.


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
MVWNtest_out <- MVWNtest(x = WN, maxlag = 10, printResults = 1)
x <- FracDiff(x = WN, d = - 0.5)
MVWNtest_out <- MVWNtest(x = x, maxlag = 10, printResults = 1)
WN_x_d <- FracDiff(x, d = 0.5)
MVWNtest_out <- MVWNtest(x = WN_x_d, maxlag = 10, printResults = 1)


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
