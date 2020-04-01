################################################################################
# 
# More Examples of FCVAR Analysis
# 
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 16, 2020
# 
################################################################################
# 
# This script demonstrates some of the additional features of the FCVAR
#   software package:
# 
#   - Forecasting
#   - Bootstrap test of hypothesis on model coefficients
#   - Bootstrap rank test
#   - Simulation of fractionally cointegrated process
# 
#  This script calls the estimates, option settings and data from
#   replication_JNP2014.m, so that file should be run before the examples
#   presented here.
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
# library(fracdiff)


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
opt$constrained  <- 0 # impose restriction dbMax ><- d ><- b ><- dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 1 # impose restriction d<-b ? 1 <- yes, 0 <- no.
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
# ESTIMATE RESTRICTED MODEL
################################################################################

DefaultOpt$gridSearch <- 0	# turn off grid search for restricted models
#	because it is too intensive.

#--------------------------------------------------------------------------------
# Test restriction that interest-rate is long-run exogenous.
#--------------------------------------------------------------------------------

k <- 2

r <- 1

opt1 <- DefaultOpt
opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)

m1r4 <- FCVARestn(x1, k, r, opt1) # This restricted model is now in the structure m1r4.


################################################################################
# FORECAST
################################################################################


# Forecast from the final restricted model.
NumPeriods <- 12 # forecast horizon set to 12 months ahead.

# Assign the model whose coefficients will be used for forecasting.
modelF <- m1r4

xf <- FCVARforecast(x1, modelF, NumPeriods)


#--------------------------------------------------------------------------------
# Testing version:
#--------------------------------------------------------------------------------
source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')

x <- x1
model <- m1r4
xf <- FCVARforecast(x, model, NumPeriods)
#--------------------------------------------------------------------------------



# Series including forecast.
seriesF <- as.matrix(rbind(x1, xf) )

# Equilibrium relation including forecasts.
equilF <- seriesF %*% modelF$coeffs$betaHat 


#--------------------------------------------------------------------------------
# Plot the series and forecast.
#--------------------------------------------------------------------------------

# Determine the size of the vertical line to delimit data and forecast
#   values.
T <- nrow(x1)
yMaxS  <- max(seriesF)
yMinS  <- min(seriesF)
yMaxEq <- max(equilF)
yMinEq <- min(equilF)

# Plot the series and forecast.
color_list <- rainbow(ncol(seriesF))
col_num <- 1
plot(seriesF[, col_num], 
     main = 'Series, including Forecast', 
     xlab = 'Time, t', 
     ylab = 'Series',
     ylim = c(yMinS, yMaxS),
     type = 'l',
     col = color_list[col_num])
abline(v = T, col = 'black', lwd = 3)
for (col_num in 2:ncol(seriesF)) {
  lines(seriesF[, col_num],
        col = color_list[col_num])
}


# Plot the equilibrium relation including forecasts.
plot(equilF, 
     main = 'Equilibrium Relation, including Forecast', 
     xlab = 'Time, t', 
     ylab = 'Equilibrium Relation',
     ylim = c(yMinEq, yMaxEq),
     type = 'l',
     col = 'black')




################################################################################
# BOOTSTRAP HYPOTHESIS TEST 
################################################################################


# Test restriction that political variables do not enter the  
#   cointegrating relation(s).

# Turn off plots for bootstrapping.
DefaultOpt$plotRoots <- 0

# Define estimation options for unrestricted model (alternative)
optUNR <- DefaultOpt

# Define estimation options for restricted model (null)
optRES <- DefaultOpt
optRES$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)

# Number of bootstrap samples to generate
B <- 999
B <- 9

# Call to open the distributed processing (comment out if unavailable)
# matlabpool ('open',4) # for versions 2013a and earlier.
# parpool # for versions 2013b and later.

source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')


# [LRbs, H, mBS, mUNR] <- FCVARboot(x1, k, r, optRES, optUNR, B)
FCVARboot_out <- FCVARboot(x1, k, r, optRES, optUNR, B)

LRbs <- FCVARboot_out$LRbs
H <- FCVARboot_out$H
mBS <- FCVARboot_out$mBS
mUNR <- FCVARboot_out$mUNR


#--------------------------------------------------------------------------------
# Compare the bootstrap distribution to chi-squared distribution
#--------------------------------------------------------------------------------

# Ispect histogram for troubleshooting.
# hist(LRbs, 
#      main = 'bootstrap density with chi-squared density')


# Estimate kernel density
LRbs_density <- density(LRbs)


# Plot bootstrap density with chi-squared density
plot(LRbs_density, 
     main = c('Bootstrap Density with Chi-squared Density', 
              sprintf('(%d bootstrap samples and %d d.f.)', 
                      B, H$df)), 
     xlab = 'Likelihood Ratio Statistic', 
     lwd = 2, 
     col = 'blue') 
LR_range <- seq(min(LRbs_density$x), max(LRbs_density$x), by = 0.01)
lines(LR_range,
      dchisq(LR_range, df = H$df), 
      col = 'red', 
      lwd = 2)
legend('topright', 
       c('Bootstrap', 'Chi-squared'), 
       # pch = 16, 
       fill = c('blue','red'))



# Estimate kernel density
# [F,XI]=ksdensity(LRbs)

# # Plot bootstrap density with chi-squared density
# figure plot(XI,F, XI, chi2pdf(XI,H.df))
# 
# legend(['Bootstrap PDF with ', num2str(B), ' BS samples'],...
#        ['Chi Squared with ', num2str(H.df),' df'])





################################################################################
# BOOTSTRAP RANK TEST 
################################################################################


# Test rank 0 against rank 1
r1 <- 0
r2 <- 1

# Number of bootstrap samples to generate
B <- 999
B <- 9

# [LR_Rnk, H_Rnk, mBSr1, mBSr2] <- FCVARbootRank(x1, k, DefaultOpt, r1, r2, B)


source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')

FCVARbootRank_out <- FCVARbootRank(x1, k, DefaultOpt, r1, r2, B)

LR_Rnk <- FCVARbootRank_out$LRbs
H_Rnk <- FCVARbootRank_out$H
mBSr1 <- FCVARbootRank_out$mBS
mBSr2 <- FCVARbootRank_out$mUNR



#--------------------------------------------------------------------------------
# Compare to P-value based on asymptotic distribution
#--------------------------------------------------------------------------------

rankTestStats <- RankTests(x1, k, DefaultOpt)

cat(sprintf('P-value: \t %1.3f\n', rankTestStats$pv[1]))

# Close distributed processing (comment out if unavailable)
# matlabpool close




################################################################################
# SIMULATION 
################################################################################


# Simulate the final restricted model, the same one used for forecasting
#   above.

# Number of periods to simulate
T_sim <- 100


# Simulate data
xSim <- FCVARsim(x1, modelF, T_sim)


#--------------------------------------------------------------------------------
# Testing version:
#--------------------------------------------------------------------------------
source('EstOptions.R')
source('FCVAR_estn.R')
source('FCVAR_lower.R')
source('FCVAR_higher.R')

x <- x1
model <- m1r4
# Simulate data
xSim <- FCVARsim(x, model, T_sim)
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# Plot the simulated series
#--------------------------------------------------------------------------------

# plot(xSim)
# legend('Support', 'Unemployment', 'Interest rate')

yMaxS  <- max(xSim)
yMinS  <- min(xSim)

# Plot the series and forecast.
color_list <- rainbow(ncol(xSim))
col_num <- 1
plot(xSim[, col_num], 
     main = 'Series, including Forecast', 
     xlab = 'Time, t', 
     ylab = 'Series',
     ylim = c(yMinS, yMaxS),
     type = 'l',
     col = color_list[col_num])
abline(v = T, col = 'black', lwd = 3)
for (col_num in 2:ncol(seriesF)) {
  lines(seriesF[, col_num],
        col = color_list[col_num])
}
legend('topleft', 
       c('Support', 'Unemployment', 'Interest Rate'), 
       # pch = 16, 
       fill = color_list)






################################################################################

endProg <- Sys.time() # stop timer
print('Total computation time from start to end:')
print(endProg - startProg) # report elapsed time

################################################################################
# End
################################################################################

