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
# Estimation
################################################################################

# FCVARoptions

opt <- FCVARoptions()

# FCVARoptionUpdates(opt, p, r)

opt <- FCVARoptions()
opt$gridSearch <- 0 # Skip grid search in optimization.
newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1)

# GetBounds(opt)

opt <- FCVARoptions()
UB_LB_bounds <- GetBounds(opt)


# FCVARestn(x,k,r,opt)

x <- votingJNP2014[, c("lib", "ir_can", "un_can")]



################################################################################
# End
################################################################################

