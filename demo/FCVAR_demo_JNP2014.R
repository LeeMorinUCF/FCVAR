################################################################################
#
# Example of Analysis Using the FCVAR Model
#
#
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business
# University of Central Florida
#
# July 28, 2021
#
################################################################################
#
# This code replicates statistics for Table 4: FCVAR results for Model 1, in
# Maggie E.C. Jones, Morten \Orregaard Nielsen & Michal Ksawery Popiel (2014).
#   "A fractionally cointegrated VAR analysis of economic voting and political support,"
#   Canadian Journal of Economics.
#
# This script also demonstrates some of the additional features of the FCVAR
#   software package:
#
#   - Forecasting
#   - Bootstrap test of hypothesis on model coefficients
#   - Bootstrap rank test
#   - Simulation of fractionally cointegrated process
#   - Grid search for starting values to improve optimization
#
# This demo script serves as the set of examples for the vignette to accompany
#   the R package FCVAR.
#
################################################################################

# Clear workspace.
rm(list = ls(all = TRUE))

# Install and attach package.
install.packages('FCVAR')
library(FCVAR)


################################################################################
# Import Data
################################################################################

x1 <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# lib is the aggregate support for the Liberal party,
# ir_can is the Canadian 3-month T-bill rate, and
# un_can is the Canadian unemployment rate.


################################################################################
# INITIALIZATION
################################################################################

p                <- ncol(x1) # system dimension.
kmax             <- 3    # Maximum number of lags for VECM.
order            <- 12   # Number of lags for white noise test in lag selection.
printWNtest      <- 1    # Print results of MVWN test to screen.

#--------------------------------------------------------------------------------
# Choosing estimation options
#--------------------------------------------------------------------------------

# Define variable to store Estimation Options (in an FCVARoptions object).
opt <- FCVARoptions(
  unrConstant = 0,
  rConstant = 0,
  levelParam = 1,
  constrained = 0,
  restrictDB = 1,
  dbMin = c(0.01, 0.01),
  dbMax = c(2.00, 2.00)
)
opt$db0 <- c(0.80, 0.80)
opt$gridSearch <- 0


################################################################################
# LAG SELECTION
################################################################################

FCVARlagSelectStats <- FCVARlagSelect(x1, kmax, p, order, opt)
# Select lag k <- 2 with minimum value of AIC.


################################################################################
# COINTEGRATION RANK TESTING
################################################################################

k <- 2
rankTestStats <- FCVARrankTests(x1, k, opt)


################################################################################
# UNRESTRICTED MODEL ESTIMATION
################################################################################

# Set parameters from specification decisions.
r <- 1
opt1 <- opt
# opt1$gridSearch <- 1
opt1$gridSearch <- 0
opt1$plotLike <- 1


# Estimate model and store in an FCVARmodel object.
m1 <- FCVARestn(x1, k, r, opt1)


# Conduct MVWN test on residuals.
MVWNtest_m1 <- MVWNtest(m1$Residuals, order, printWNtest)


################################################################################
# IMPOSE RESTRICTIONS AND TEST THEM
################################################################################


#--------------------------------------------------------------------------------
# Test restriction that d = b = 1.
#--------------------------------------------------------------------------------

# Set options to impose restriction.
opt1 <- opt
opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt1$r_psi <- 1


# Estimate model and store in an FCVARmodel object.
m1r1 <- FCVARestn(x1, k, r, opt1)


# Conduct MVWN test on residuals.
MVWNtest_m1r1 <- MVWNtest(m1r1$Residuals, order, printWNtest)


# Test the null of m1r1 against the alternative m1.
Hdb <- FCVARhypoTest(m1, m1r1)


#--------------------------------------------------------------------------------
# Test restriction that political variables do not enter the cointegrating relation(s).
#--------------------------------------------------------------------------------

# Set options to impose restriction.
opt1 <- opt
opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)


# Estimate model and store in an FCVARmodel object.
m1r2 <- FCVARestn(x1, k, r, opt1)


# Conduct MVWN test on residuals.
MVWNtest_m1r2 <- MVWNtest(m1r2$Residuals, order, printWNtest)


# Test the null of m1r2 against the alternative m1.
Hbeta1 <- FCVARhypoTest(m1, m1r2)


#--------------------------------------------------------------------------------
# Test restriction that political variable is long-run exogenous.
#--------------------------------------------------------------------------------

# Set options to impose restriction.
opt1 <- opt
opt1$R_Alpha <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)


# Estimate model and store in an FCVARmodel object.
m1r3 <- FCVARestn(x1, k, r, opt1)


# Conduct MVWN test on residuals.
MVWNtest_m1r3 <- MVWNtest(m1r3$Residuals, order, printWNtest)


# Test the null of m1r3 against the alternative m1
Halpha1 <- FCVARhypoTest(m1, m1r3)


#--------------------------------------------------------------------------------
# Test restriction that interest-rate is long-run exogenous.
#--------------------------------------------------------------------------------

# Set options to impose restriction.
opt1 <- opt
opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)


# Estimate model and store in an FCVARmodel object.
m1r4 <- FCVARestn(x1, k, r, opt1)


# Conduct MVWN test on residuals.
MVWNtest_m1r4 <- MVWNtest(m1r4$Residuals, order, printWNtest)


# Test the null of m1r4 against the alternative m1.
Halpha2 <- FCVARhypoTest(m1, m1r4)


#--------------------------------------------------------------------------------
# Test restriction that unemployment is long-run exogenous.
#--------------------------------------------------------------------------------

# Set options to impose restriction.
opt1 <- opt
opt1$R_Alpha <- matrix(c(0, 0, 1), nrow = 1, ncol = 3)


# Estimate model and store in an FCVARmodel object.
m1r5 <- FCVARestn(x1, k, r, opt1)


# Conduct MVWN test on residuals.
MVWNtest_m1r5 <- MVWNtest(m1r5$Residuals, order, printWNtest)


# Test the null of m1r5 against the alternative m1.
Halpha3 <- FCVARhypoTest(m1, m1r5)


#--------------------------------------------------------------------------------
# RESTRICTED MODEL OUTPUT
#   - print normalized beta and alpha for model m1r4.
#--------------------------------------------------------------------------------

# Assign model.
modelRstrct <- m1r4

# Perform Normalization.
G <- solve(modelRstrct$coeffs$betaHat[1:r, 1:r])
betaHatR <- modelRstrct$coeffs$betaHat %*% G
# alphaHat is post multiplied by G^{-1} so that pi = a(G^{-1})Gb' = ab'
alphaHatR <- modelRstrct$coeffs$alphaHat %*% t(solve(G))


# Print output.
print("betaHatR' = ")
print(t(betaHatR), print.gap = 5)
print("alphaHatR' = ")
print(t(alphaHatR), print.gap = 5)


################################################################################
# End
################################################################################




