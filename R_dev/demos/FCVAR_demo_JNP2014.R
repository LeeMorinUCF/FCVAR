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
# opt1$gridSearch <- 0
# opt1$gridSearch <- 1
# opt1$plotLike <- 1


# Estimate model and store in an FCVARmodel object.
m1 <- FCVARestn(x1, k, r, opt1)


# Conduct MVWN test on residuals.
MVWNtest_m1 <- MVWNtest(m1$Residuals, order, printWNtest)


# Experiment by imposing linear restriction d = b
# through restriction matrices.
opt1 <- opt
opt1$restrictDB <- 0
# Alternate approach to impose d = b:
opt1$R_psi <- matrix(c(1, -1), nrow = 1, ncol = 2)
# Other example with linear restriction on d and b:
# opt1$R_psi <- matrix(c(1, -2), nrow = 1, ncol = 2)
opt1$r_psi <- 0
m1 <- FCVARestn(x1, k, r, opt1)
# Error in R %*% Hinv : non-conformable arguments
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# R> warnings()
# Warning messages:
# 1: In cbind(param, matrix(coeffs$muHat, nrow = 1, ncol = p)) :
#   number of rows of result is not a multiple of vector length (arg 1)
# 2: In matrix(param[1:length(param)], nrow = p, ncol = p *  ... :
#     data length [17] is not a sub-multiple or multiple of the number of rows [3]
# 3: In matrix(param[1:length(param)], nrow = p, ncol = p *  ... :
#     data length [17] is not a sub-multiple or multiple of the number of rows [3]

# Solution (a change within function SEmat2vecU in script FCVAR_aux.R):
# # Append the level parameter MuHat.
# if (opt$levelParam) {
#
#   # print('Append the level parameter MuHat.')
#   # print('opt$levelParam = ')
#   # print(opt$levelParam)
#   # print('param = ')
#   # print(param)
#   # print('coeffs$muHat = ')
#   # print(coeffs$muHat)
#
#   # param <- cbind(param, matrix(coeffs$muHat, nrow = 1, ncol = p))
#   param <- cbind(matrix(param, nrow = 1, ncol = length(param)),
#                  matrix(coeffs$muHat, nrow = 1, ncol = p))
#
#   # print('new_param = ')
#   # print(param)
# }
#
# And similarly for the others.


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
# Forecasting with the FCVAR Model
################################################################################


# Forecast from the final restricted model.
NumPeriods <- 12 # forecast horizon set to 12 months ahead.

# Assign the model whose coefficients will be used for forecasting.
modelF <- m1r4

# Calculate the forecasts.
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
T_xf <- nrow(x1)
yMaxS  <- max(seriesF)
yMinS  <- min(seriesF)
yMaxEq <- max(equilF)
yMinEq <- min(equilF)

# Plot the series and forecast.
# fig_ext <- 'png'
# out_file_path <- sprintf('%s/forecast_vars.%s', fig_dir, fig_ext)
color_list <- rainbow(ncol(seriesF))
label_list <- c('Liberal Support', 'CDN 3-mo T-Bill', 'CDN Unemp. Rate')
lty_list <- c('solid', 'dashed', 'dot')
col_num <- 1
# png(out_file_path)
plot(seriesF[, col_num],
     main = 'Series, including Forecast',
     xlab = 'Time, t',
     ylab = 'Series',
     ylim = c(yMinS, yMaxS),
     type = 'l',
     lwd = 3,
     col = color_list[col_num])
abline(v = T_xf, col = 'black', lwd = 3, lty = 'dashed')
for (col_num in 2:ncol(seriesF)) {
  lines(seriesF[, col_num],
        lwd = 3,
        lty = col_num,
        col = color_list[col_num])
}
# Settings in GUI:
legend(197, 22.75, legend = label_list,
       col = color_list, lty = 1:3, cex = 0.65)
# Settings in article:
# legend(150, 20, legend = label_list,
#        col = color_list, lty = 1:3, cex = 1.0)
# dev.off()


# Plot the equilibrium relation including forecasts.
# out_file_path <- sprintf('%s/forecast_eqbm.%s', fig_dir, fig_ext)
# png(out_file_path)
plot(equilF,
     main = 'Equilibrium Relation, including Forecast',
     xlab = 'Time, t',
     ylab = 'Equilibrium Relation',
     ylim = c(yMinEq, yMaxEq),
     type = 'l',
     lwd = 3,
     col = 'black')
abline(v = T_xf, col = 'black', lwd = 3, lty = 'dashed')
# dev.off()


#--------------------------------------------------------------------------------
# Plot the series and forecast together.
#--------------------------------------------------------------------------------

# orig_par <- par()

# out_file_path <- sprintf('%s/forecast_vars_eqbm.%s', fig_dir, fig_ext)
plot.new()
# png(out_file_path)

# Bottom pane: Equilibrium Relation
par(fig = c(0, 1, 0, 0.4), new = TRUE, mar = c(5.1, 4.1, 1.1, 2.1))
# Reset plotting parameters after.
# par(fig = c(0, 1, 0, 1), new = FALSE, mar = c(5.1, 4.1, 4.1, 2.1))
plot(equilF,
     # main = 'Equilibrium Relation, including Forecast',
     xlab = 'Time, t',
     # ylab = 'Equilibrium Relation',
     ylab = 'Equilibrium',
     ylim = c(yMinEq, yMaxEq),
     type = 'l',
     lwd = 3,
     col = 'black')
abline(v = T_xf, col = 'black', lwd = 3, lty = 'dashed')

# Top pane: Series and Forecast.
par(fig = c(0, 1, 0.4, 1.0), new = TRUE, mar = c(2.1, 4.1, 4.1, 2.1))
col_num <- 1
plot(seriesF[, col_num],
     # main = 'Series, including Forecast',
     main = 'Variables and Equilibrium Relation with Forecast',
     # xlab = 'Time, t',
     ylab = 'Variables',
     ylim = c(yMinS, yMaxS),
     type = 'l',
     lwd = 3,
     col = color_list[col_num])
abline(v = T_xf, col = 'black', lwd = 3, lty = 'dashed')
for (col_num in 2:ncol(seriesF)) {
  lines(seriesF[, col_num],
        lwd = 3,
        lty = col_num,
        col = color_list[col_num])
}
# Settings in article:
# legend(200, 20, legend = label_list,
#        col = color_list, lty = 1:3, cex = 1.0)
# Settings in GUI:
legend(220, 22, legend = label_list,
       col = color_list, lty = 1:3, cex = 0.55)
# dev.off()

# Reset plotting parameters after plot is complete.
par(fig = c(0, 1, 0, 1), new = FALSE, mar = c(5.1, 4.1, 4.1, 2.1))


################################################################################
# Bootstrap Hypothesis Test
################################################################################


# Test restriction that political variables do not enter the
#   cointegrating relation(s).

# Turn off plots for bootstrapping.
opt$plotRoots <- 0

# Define estimation options for unrestricted model (alternative)
optUNR <- opt

# Define estimation options for restricted model (null)
optRES <- opt
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
# out_file_name <- 'R_dev/FCVARboot.RData'
# save.image(file = out_file_name)

#--------------------------------------------------------------------------------
# Compare the bootstrap distribution to chi-squared distribution
#--------------------------------------------------------------------------------

# Estimate kernel density
LRbs_density <- density(LRbs)


# Plot bootstrap density with chi-squared density.
# out_file_path <- sprintf('%s/LRdensity_bw045.%s', fig_dir, fig_ext)
# png(out_file_path)
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
# dev.off()



################################################################################
# Bootstrap Rank Test
################################################################################


# Test rank 0 against rank 1
r1 <- 0
r2 <- 1

# Number of bootstrap samples to generate
B <- 999

set.seed(42)
FCVARbootRank_stats <- FCVARbootRank(x1, k, opt, r1, r2, B)

LR_Rnk <- FCVARbootRank_stats$LRbs
H_Rnk <- FCVARbootRank_stats$H
mBSr1 <- FCVARbootRank_stats$mBS
mBSr2 <- FCVARbootRank_stats$mUNR


# # Save results.
# out_file_name <- 'R_dev/FCVARbootRank.RData'
# save.image(file = out_file_name)


#--------------------------------------------------------------------------------
# Compare to P-value based on asymptotic distribution
#--------------------------------------------------------------------------------

rankTestStats <- FCVARrankTests(x1, k, opt)

cat(sprintf('P-value (asy): \t %1.3f\n', rankTestStats$pv[1]))






################################################################################
# Simulating from the FCVAR Model
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
# out_file_path <- sprintf('%s/sim.%s', fig_dir, fig_ext)
# png(out_file_path)
color_list <- rainbow(ncol(x_sim))
col_num <- 1
plot(x_sim[, col_num],
     main = 'Simulated Data',
     xlab = 'Time, t',
     ylab = 'Series',
     ylim = c(yMinS, yMaxS),
     type = 'l',
     lwd = 3,
     col = color_list[col_num])
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
legend('topleft', legend = label_list,
       col = color_list, lty = 1:3, cex = 0.65)
# Settings in article:
# legend('topleft', legend = label_list,
#        col = color_list, lty = 1:3, cex = 1.0)
# dev.off()




################################################################################
# Grid Search for Local Optima
################################################################################

# Restrict equality of fractional parameters.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.01
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 1 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$progress     <- 2 # Show progress report on each value of b.
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x1, k = 2, r = 1, opt)

# Output plot for article.
# out_file_path <- sprintf('%s/gridDB.%s', fig_dir, fig_ext)
# plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt,
#                     main = 'default',
#                     file = out_file_path,
#                     file_ext = fig_ext)
# out_file_path <- 'R_dev/Figures/gridDB.png'
# grDevices::png(out_file_path)
# graphics::plot(x = likeGrid_params)
# grDevices::dev.off()


# Linear restriction on fractional parameters.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.01
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
# Impose linear restriction on d and b:
opt$R_psi        <- matrix(c(2, -1), nrow = 1, ncol = 2)
opt$r_psi        <- 0.5
opt$progress     <- 2 # Show progress report on each value of b.
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x1, k = 2, r = 1, opt)

# Output plot for article.
# out_file_path <- sprintf('%s/gridPhi.%s', fig_dir, fig_ext)
# plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt,
#                     main = 'default',
#                     file = out_file_path,
#                     file_ext = fig_ext)
# out_file_path <- 'R_dev/Figures/gridPhi.png'
# grDevices::png(out_file_path)
# graphics::plot(x = likeGrid_params)
# grDevices::dev.off()


# Constrained 2-dimensional optimization.
# Impose restriction dbMax >= d >= b >= dbMin.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.01
opt$dbStep2D     <- 0.02
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 1 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$progress     <- 2 # Show progress report on each value of b.
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x1, k = 2, r = 1, opt)

# Output plot for article.
# out_file_path <- sprintf('%s/gridConst.%s', fig_dir, fig_ext)
# plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt,
#                     main = 'default',
#                     file = out_file_path,
#                     file_ext = fig_ext)
# out_file_path <- 'R_dev/Figures/gridConst.png'
# grDevices::png(out_file_path)
# graphics::plot(x = likeGrid_params)
# grDevices::dev.off()


# Unconstrained 2-dimensional optimization.
opt <- FCVARoptions()
opt$dbStep1D     <- 0.01
opt$dbStep2D     <- 0.01
# opt$dbStep2D     <- 0.2 # Faster with a coarser grid.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
opt$progress     <- 2 # Show progress report on each value of b.
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
likeGrid_params <- FCVARlikeGrid(x1, k = 2, r = 1, opt)

# Output plot for article.
# out_file_path <- sprintf('%s/grid3d.%s', fig_dir, fig_ext)
# plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt,
#                     main = 'default',
#                     file = out_file_path,
#                     file_ext = fig_ext)
# out_file_path <- 'R_dev/Figures/grid3d.png'
# grDevices::png(out_file_path)
# graphics::plot(x = likeGrid_params)
# grDevices::dev.off()






################################################################################
# End
################################################################################




