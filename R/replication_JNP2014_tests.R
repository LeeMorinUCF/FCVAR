################################################################################
# 
# Example of FCVAR analysis
# 
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 10, 2020
# 
################################################################################
# 
# This code replicates Table 4: FCVAR results for Model 1, in
# Maggie E.C. Jones, Morten \Orregaard Nielsen & Michal Ksawery Popiel (2014).
#   "A fractionally cointegrated VAR analysis of economic voting and political support,"
#   Canadian Journal of Economics.
# 
################################################################################


################################################################################
# Clearing Workspace and Declaring Packages
################################################################################

# Clear workspace.
rm(list = ls(all = TRUE))


################################################################################
# Set parameters for file IO
################################################################################

# Set working directory.
wd_path <- '~/Research/FCVAR/GitRepo/FCVAR/R'
setwd(wd_path)

# Add path containing Auxillary files required for estimation.
# aux_path <- '~/Research/FCVAR/GitRepo/FCVAR/R/Auxiliary'
# Load library containing Auxillary functions required for estimation.
# source('Auxiliary.R')

# Load library containing Estimation Options required for estimation.
source('EstOptions.R')
# Contains a single function EstOptions() with default options. 

# Load library with Main Estimation Function for FCVAR
source('FCVAR_estn.R')

# Load library of Functions for FCVAR Estimation
source('FCVAR_lower.R')
# Load library of Pre- and Post-Estimation Functions for FCVAR
source('FCVAR_higher.R')

# Load fracdiff library for testing.
library(fracdiff)


################################################################################
# Import Data
################################################################################


data <- read.csv('data_JNP2014.csv')

# data for each model.
x1 <- data[, c(1, 3, 5)]
x2 <- data[, c(2, 3, 5)]
x3 <- data[, c(1, 2, 3, 5)]
x4 <- data[, c(1, 3, 4, 5, 6)]
x5 <- data[, c(2, 3, 4, 5, 6)]
x6 <- data[, c(1, 2, 3, 4, 5, 6)]


################################################################################
# INITIALIZATION
################################################################################

p                <- ncol(x1) # system dimension.
kmax             <- 3    # maximum number of lags for VECM.
order            <- 12   # number of lags for white noise test in lag selection.
printWNtest      <- 1    # to print results of white noise tests post-estimation.

#--------------------------------------------------------------------------------
# Choosing estimation options
#--------------------------------------------------------------------------------

opt <- EstOptions() # Define variable to store Estimation Options (object).

opt$dbMin        <- c(0.01, 0.01) # lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # upper bound for d,b.
opt$unrConstant  <- 0 # include an unrestricted constant? 1 <- yes, 0 <- no.
opt$rConstant    <- 0 # include a restricted constant? 1 <- yes, 0 <- no.
opt$levelParam   <- 1 # include level parameter? 1 <- yes, 0 <- no.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 1 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$db0          <- c(0.8, 0.8) # set starting values for optimization algorithm.
opt$N            <- 0 # number of initial values to condition upon.
opt$print2screen <- 1 # print output.
opt$printRoots   <- 1 # do not print roots of characteristic polynomial.
opt$plotRoots    <- 1 # do not plot roots of characteristic polynomial.
opt$gridSearch   <- 1 # For more accurate estimation, perform the grid search.
					  # This will make estimation take longer.
opt$plotLike     <- 0 # Plot the likelihood (if gridSearch <- 1).
opt$progress 	   <- 0 # Show grid search progress indicator waitbar.
opt$updateTime   <- 0.5 # How often progress is updated (seconds).

# Linux example:
opt$progLoc <- '"/usr/bin/fdpval"'  # location path with program name
                                    # of fracdist program, if installed
                                    # Note: use both single (outside) and double
                                    # quotes (inside). This is especially important
                                    # if path name has spaces.
# Windows example:
# opt$progLoc <- '".\fdpval\fdpval"'  # program located in folder fdpval in current directory

# There are many other options (see EstOptions.m for
# everything else. These can be accessed/modified as in, for example:
# opt$dbFminOptions$Algorithm <- 'interior-point'

DefaultOpt <- opt # Store the options for restoring them in between hypothesis tests.

# startProg <- tic() # start timer
startProg <- Sys.time() # start timer


################################################################################
# LAG SELECTION 
################################################################################

# LagSelect(x1, kmax, p, order, opt)

x <- x1

source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')
opt <- DefaultOpt
opt$gridSearch <- 0
opt$db0 <- c(0.676, 0.676)
# opt$gridSearch <- 1
# opt$progress <- 2
LagSelect(x, kmax, p, order, opt)

#--------------------------------------------------------------------------------
# Testing LkeGrid()
#--------------------------------------------------------------------------------


source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')
opt <- DefaultOpt
opt$gridSearch <- 1
opt$progress <- 2
k <- 3
r <- 3
# Shouldn't this be done here, too?
# Update options based on initial user input.
opt <- updateRestrictions(opt, p, r)
LikeGrid(x, k, r, opt)


################################################################################
# COINTEGRATION RANK TESTING
################################################################################


k <- 2

# rankTestStats <- RankTests(x1, k, opt)


source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')
k <- 2
opt <- DefaultOpt
opt$gridSearch <- 0
rankTestStats <- RankTests(x, k, opt)


# Testing the gridsearch option.
source('FCVAR_lower.R')
source('FCVAR_higher.R')
opt2 <- DefaultOpt
# opt2$gridSearch <- 0
opt2$gridSearch <- 1
startRankTest <- Sys.time() # start timer

rankTestStats <- RankTests(x, k, opt2)

endRankTest <- Sys.time() # stop timer
print(endRankTest - startRankTest) # report elapsed time


################################################################################
# UNRESTRICTED MODEL ESTIMATION
################################################################################


k <- 2

r <- 1

opt1 <- DefaultOpt
opt1$gridSearch <- 0

#--------------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------------

opt <- opt1
opt$gridSearch <- 0
x <- x1

source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')

# m1 <- FCVARestn(x1, k, r, opt1) # This model is now in the structure m1.
m1 <- FCVARestn(x, k, r, opt) # This model is now in the structure m1.



source('FCVAR_higher.R')
mv_wntest_m1 <- mv_wntest(m1$Residuals, order, printWNtest)

################################################################################
# IMPOSE RESTRICTIONS AND TEST THEM
################################################################################




DefaultOpt$gridSearch <- 0	# turn off grid search for restricted models
							#	because it is too intensive.

#--------------------------------------------------------------------------------
# Test restriction that d<-b<-1.
#--------------------------------------------------------------------------------


opt1 <- DefaultOpt
opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt1$r_psi <- 1

m1r1 <- FCVARestn(x1, k, r, opt1) # This restricted model is now in the structure m1r1.


#--------------------------------------------------------------------------------
# Testing version:
#--------------------------------------------------------------------------------
source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')
opt <- DefaultOpt
opt$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt$r_psi <- 1
m1r1 <- FCVARestn(x, k, r, opt) # This restricted model is now in the structure m1r1.
#--------------------------------------------------------------------------------



mv_wntest_m1r1 <- mv_wntest(m1r1$Residuals, order, printWNtest)

Hdb <- HypoTest(m1, m1r1) 	# Test the null of m1r1 against the alternative m1 and
							# store the results in the structure Hdb.





#--------------------------------------------------------------------------------
# Test restriction that political variables do not enter the cointegrating relation(s).
#--------------------------------------------------------------------------------


opt1 <- DefaultOpt
opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)

m1r2 <- FCVARestn(x1, k, r, opt1) # This restricted model is now in the structure m1r2.


#--------------------------------------------------------------------------------
# Testing version:
#--------------------------------------------------------------------------------
source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')
opt <- DefaultOpt
opt$gridSearch <- 0
opt$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
m1r2 <- FCVARestn(x, k, r, opt) # This restricted model is now in the structure m1r2.
#--------------------------------------------------------------------------------


mv_wntest_m1r2 <- mv_wntest(m1r2$Residuals, order, printWNtest)

Hbeta1 <- HypoTest(m1, m1r2) 	# Test the null of m1r2 against the alternative m1 and
								# store the results in the structure Hbeta1.






#--------------------------------------------------------------------------------
# Test restriction that political variable is long-run exogenous.
#--------------------------------------------------------------------------------

opt1 <- DefaultOpt
opt1$R_Alpha <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
opt1$gridSearch <- 0

m1r3 <- FCVARestn(x1, k, r, opt1) # This restricted model is now in the structure m1r3.

#--------------------------------------------------------------------------------
# Testing version:
#--------------------------------------------------------------------------------
source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')
opt <- DefaultOpt
opt$R_Alpha <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
opt$gridSearch <- 0
m1r3 <- FCVARestn(x, k, r, opt) # This restricted model is now in the structure m1r3.
#--------------------------------------------------------------------------------


mv_wntest_m1r3 <- mv_wntest(m1r3$Residuals, order, printWNtest)

Halpha1 <- HypoTest(m1, m1r3) 	# Test the null of m1r3 against the alternative m1 and
								# store the results in the structure Halpha1.




#--------------------------------------------------------------------------------
# Test restriction that interest-rate is long-run exogenous.
#--------------------------------------------------------------------------------


## Test restriction that interest-rate is long-run exogenous.
opt1 <- DefaultOpt
opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
opt1$gridSearch <- 0

m1r4 <- FCVARestn(x1, k, r, opt1) # This restricted model is now in the structure m1r4.

mv_wntest_m1r4 <- mv_wntest(m1r4$Residuals, order, printWNtest)

Halpha2 <- HypoTest(m1, m1r4) 	# Test the null of m1r4 against the alternative m1 and
								# store the results in the structure Halpha2.





#--------------------------------------------------------------------------------
## Test restriction that unemployment is long-run exogenous.
#--------------------------------------------------------------------------------


opt1 <- DefaultOpt
opt1$gridSearch <- 0
opt1$R_Alpha <- matrix(c(0, 0, 1), nrow = 1, ncol = 3)
k<-2 
r <-1
m1r5 <- FCVARestn(x1, k, r, opt1) # This restricted model is now in the structure m1r5.

mv_wntest_m1r5 <- mv_wntest(m1r5$Residuals, order, printWNtest)

Halpha3 <- HypoTest(m1, m1r5) 	# Test the null of m1r5 against the alternative m1 and
								# store the results in the structure Halpha3.



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

# Yes, I know we shouldn't be inverting, Harry, but this is quick and easy.

# Print output.
print(betaHatR)
print(alphaHatR)




# toc(startProg)
endProg <- Sys.time() # stop timer
print('Total computation time from start to end:')
print(endProg - startProg) # report elapsed time



################################################################################
# End
################################################################################

