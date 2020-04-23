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
#   FCVAR_replication_JNP2014.R, so that file should be run before the examples
#   presented here.
#
################################################################################


################################################################################
# Import Data
################################################################################

#
# data <- read.csv('data_JNP2014.csv')
#
# # data for each model.
# x1 <- data[, c(1, 3, 5)]
# x2 <- data[, c(2, 3, 5)]
# x3 <- data[, c(1, 2, 3, 5)]
# x4 <- data[, c(1, 3, 4, 5, 6)]
# x5 <- data[, c(2, 3, 4, 5, 6)]
# x6 <- data[, c(1, 2, 3, 4, 5, 6)]


################################################################################
# INITIALIZATION
################################################################################

# p                <- ncol(x1) # system dimension.
# kmax             <- 3    # maximum number of lags for VECM.
# order            <- 12   # number of lags for white noise test in lag selection.
# printWNtest      <- 1    # to print results of white noise tests post-estimation.

#--------------------------------------------------------------------------------
# Choosing estimation options
#--------------------------------------------------------------------------------

# opt <- EstOptions() # Define variable to store Estimation Options (object).
#
# opt$dbMin        <- c(0.01, 0.01) # lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # upper bound for d,b.
# opt$unrConstant  <- 0 # include an unrestricted constant? 1 <- yes, 0 <- no.
# opt$rConstant    <- 0 # include a restricted constant? 1 <- yes, 0 <- no.
# opt$levelParam   <- 1 # include level parameter? 1 <- yes, 0 <- no.
# opt$constrained  <- 0 # impose restriction dbMax ><- d ><- b ><- dbMin ? 1 <- yes, 0 <- no.
# opt$restrictDB   <- 1 # impose restriction d<-b ? 1 <- yes, 0 <- no.
# opt$db0          <- c(0.8, 0.8) # set starting values for optimization algorithm.
# opt$N            <- 0 # number of initial values to condition upon.
# opt$print2screen <- 1 # print output.
# opt$printRoots   <- 1 # do not print roots of characteristic polynomial.
# opt$plotRoots    <- 1 # do not plot roots of characteristic polynomial.
# opt$gridSearch   <- 1 # For more accurate estimation, perform the grid search.
# # This will make estimation take longer.
# opt$plotLike     <- 0 # Plot the likelihood (if gridSearch <- 1).
# opt$progress 	   <- 0 # Show grid search progress indicator waitbar.
# opt$updateTime   <- 0.5 # How often progress is updated (seconds).
#
# # Linux example:
# opt$progLoc <- '"/usr/bin/fdpval"'  # location path with program name
# # of fracdist program, if installed
# # Note: use both single (outside) and double
# # quotes (inside). This is especially important
# # if path name has spaces.
# # Windows example:
# # opt$progLoc <- '".\fdpval\fdpval"'  # program located in folder fdpval in current directory
#
# # There are many other options (see EstOptions.m for
# # everything else. These can be accessed/modified as in, for example:
# # opt$dbFminOptions$Algorithm <- 'interior-point'
#
# DefaultOpt <- opt # Store the options for restoring them in between hypothesis tests.
#
# # startProg <- tic() # start timer
# startProg <- Sys.time() # start timer

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

m1r4 <- FCVARestn(x1, k, r, opt1)


################################################################################
# FORECAST
################################################################################


# Forecast from the final restricted model.
NumPeriods <- 12 # forecast horizon set to 12 months ahead.

# Assign the model whose coefficients will be used for forecasting.
modelF <- m1r4

xf <- FCVARforecast(x1, modelF, NumPeriods)


# Append forecast to series.
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
fig_ext <- 'png'
out_file_path <- sprintf('R_dev/Figures/forecast_vars.%s', fig_ext)
color_list <- rainbow(ncol(seriesF))
label_list <- c('Liberal Support', 'CDN 3-mo T-Bill', 'CDN Unemp. Rate')
lty_list <- c('solid', 'dashed', 'dot')
col_num <- 1
png(out_file_path)
plot(seriesF[, col_num],
     main = 'Series, including Forecast',
     xlab = 'Time, t',
     ylab = 'Series',
     ylim = c(yMinS, yMaxS),
     type = 'l',
     lwd = 3,
     col = color_list[col_num])
abline(v = T, col = 'black', lwd = 3, lty = 'dashed')
for (col_num in 2:ncol(seriesF)) {
  lines(seriesF[, col_num],
        lwd = 3,
        lty = col_num,
        col = color_list[col_num])
}
# Settings in GUI:
# legend(197, 22.75, legend = label_list,
#        col = color_list, lty = 1:3, cex = 0.65)
# Settings in article:
legend(150, 20, legend = label_list,
       col = color_list, lty = 1:3, cex = 1.0)
dev.off()


# Plot the equilibrium relation including forecasts.
out_file_path <- sprintf('R_dev/Figures/forecast_eqbm.%s', fig_ext)
png(out_file_path)
plot(equilF,
     main = 'Equilibrium Relation, including Forecast',
     xlab = 'Time, t',
     ylab = 'Equilibrium Relation',
     ylim = c(yMinEq, yMaxEq),
     type = 'l',
     lwd = 3,
     col = 'black')
abline(v = T, col = 'black', lwd = 3, lty = 'dashed')
dev.off()




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
# B <- 99

set.seed(42)
FCVARboot_stats <- FCVARboot(x1, k, r, optRES, optUNR, B)

LRbs <- FCVARboot_stats$LRbs
H <- FCVARboot_stats$H
mBS <- FCVARboot_stats$mBS
mUNR <- FCVARboot_stats$mUNR


# Save results.
out_file_name <- 'R_dev/FCVARboot.RData'
save.image(file = out_file_name)

#--------------------------------------------------------------------------------
# Compare the bootstrap distribution to chi-squared distribution
#--------------------------------------------------------------------------------

# Estimate kernel density
LRbs_density <- density(LRbs)


# Plot bootstrap density with chi-squared density.
out_file_path <- sprintf('R_dev/Figures/LRdensity_bw045.%s', fig_ext)
png(out_file_path)
plot(LRbs_density,
     main = c('Bootstrap Density with Chi-squared Density',
              sprintf('(%d bootstrap samples and %d d.f.)',
                      B, H$df)),
     xlab = 'Likelihood Ratio Statistic',
     xlim = c(-20, 25),
     ylim = c(0, 0.5),
     lwd = 3,
     col = 'blue')
LR_range <- seq(min(LRbs_density$x), max(LRbs_density$x), by = 0.01)
lines(LR_range,
      dchisq(LR_range, df = H$df),
      col = 'red',
      lwd = 3,
      lty = 'dotted')
legend('topleft',
       c('Bootstrap', 'Chi-squared'),
       # pch = 16,
       col = c('blue','red'),
       # fill = c('blue','red'),
       cex = 1.0,
       lty = c(1, 3))
dev.off()




################################################################################
# BOOTSTRAP RANK TEST
################################################################################


# Test rank 0 against rank 1
r1 <- 0
r2 <- 1

# Number of bootstrap samples to generate
B <- 999
# B <- 9

set.seed(42)
FCVARbootRank_stats <- FCVARbootRank(x1, k, DefaultOpt, r1, r2, B)

LR_Rnk <- FCVARbootRank_stats$LRbs
H_Rnk <- FCVARbootRank_stats$H
mBSr1 <- FCVARbootRank_stats$mBS
mBSr2 <- FCVARbootRank_stats$mUNR


# Save results.
out_file_name <- 'R_dev/FCVARbootRank.RData'
save.image(file = out_file_name)


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
set.seed(42)
x_sim <- FCVARsim(x1, modelF, T_sim)



#--------------------------------------------------------------------------------
# Plot the simulated series
#--------------------------------------------------------------------------------

yMaxS  <- max(x_sim)
yMinS  <- min(x_sim)

# Plot the series and forecast.
out_file_path <- sprintf('R_dev/Figures/sim.%s', fig_ext)
png(out_file_path)
color_list <- rainbow(ncol(xSim))
col_num <- 1
plot(x_sim[, col_num],
     main = 'Simulated Data',
     xlab = 'Time, t',
     ylab = 'Series',
     ylim = c(yMinS, yMaxS),
     type = 'l',
     lwd = 3,
     col = color_list[col_num])
# abline(v = T, col = 'black', lwd = 3)
for (col_num in 2:ncol(x_sim)) {
  lines(x_sim[, col_num],
        lwd = 3,
        lty = col_num,
        col = color_list[col_num])
}
# Original:
# legend('topleft',
#        c('Support', 'Unemployment', 'Interest Rate'),
#        # pch = 16,
#        fill = color_list)
# Settings in GUI:
# legend(197, 22.75, legend = label_list,
#        col = color_list, lty = 1:3, cex = 0.65)
# Settings in article:
legend('topleft', legend = label_list,
       col = color_list, lty = 1:3, cex = 1.0)
dev.off()




################################################################################
# End
################################################################################

