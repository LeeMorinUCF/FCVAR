################################################################################
#
# Examples of Underidentified Restrictions in the FCVAR Model
#
#
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business
# University of Central Florida
#
# May 25, 2021
#
################################################################################
#
# This is a series of examples of tests of restrictions in the FCVAR package.
# Some of the examples throw errors because the restricted model is underidentified.
#
################################################################################

################################################################################
# Load data for example.
################################################################################


fcvar_path <- "C:/Users/le279259/OneDrive - University of Central Florida/Documents/Research/FCVAR"
ex_path <- sprintf("%s/Restriction_Example/Brent_Prices.csv", fcvar_path)
brent <- read.csv(ex_path)

summary(brent)

plot(brent[, 'SPOT'])

# Select subset of data.
x1 <- brent[,seq(2, 4)]

summary(x1)




################################################################################
# Pre-estimation for model specification
################################################################################

opt <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(1.20, 1.20), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  restrictDB   = 0
)

p                <- ncol(x1) # system dimension.
kmax             <- 4    # maximum number of lags for VECM.
order            <- 12   # number of lags for white noise test in lag selection.



FCVARlagSelectStats <- FCVARlagSelect(x = x1, kmax, r = 3, order, opt)


k <- 2

rankTestStats <- FCVARrankTests(x = x1, k, opt)




################################################################################
# Estimate standard models.
################################################################################


m1 <- FCVARestn(x = x1, k = 2, r = 2, opt)
# Error in solve.default(H) :
#   Lapack routine dgesv: system is exactly singular: U[5,5] = 0
# 3.
# solve.default(H)
# 2.
# solve(H) at FCVAR_estn.R#552
# 1.
# FCVARestn(x = x1, k = 2, r = 2, opt)



m1 <- FCVARestn(x = x1, k = 2, r = 1, opt)






################################################################################
# Estimate underidentified models.
################################################################################





################################################################################
# End
################################################################################
