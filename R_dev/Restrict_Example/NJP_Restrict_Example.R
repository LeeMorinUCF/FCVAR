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

opt <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  plotRoots    = 0 # Don't create plots for tests.
)
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
results <- FCVARestn(x, k = 2, r = 1, opt)



opt1 <- opt
opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt1$r_psi <- 1
m1r1 <- FCVARestn(x, k = 2, r = 1, opt1)



# Relax the constraint d = b while restricting d = 1.

# With or without the d = 1 restriction, there is a warning
# about incompatible dimensions.

# Be careful with this one.
opt1a <- opt
opt1a$restrictDB <- 0
opt1a$constrained <- 0
opt1a$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
opt1a$r_psi <- 1
m1r1a <- FCVARestn(x, k = 2, r = 1, opt1a)



# A restricted model.


# Or start a new object with otherwise same restrictions.
opt1 <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  R_Beta       = matrix(c(1, 0, 0), nrow = 1, ncol = 3),
  plotRoots    = 0 # Don't create plots for tests.
)
m1r2 <- FCVARestn(x, k = 2, r = 1, opt1)


# Add a degree of freedom.

opt <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  plotRoots    = 0 # Don't create plots for tests.
)
m2 <- FCVARestn(x, k = 2, r = 2, opt)



# Test equality restrictions on beta with rank 2.
opt2 <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  # R_Beta       = matrix(c(1, 0, 0), nrow = 1, ncol = 3),
  # R_Beta       = cbind(diag(4),matrix(0, nrow=4, ncol=6)),
  # R_Beta       = cbind(diag(4),matrix(0, nrow=4, ncol=2)),
  R_Beta       = matrix(c(1, 0, 0, 0, 0, 0,
                          0, 1, 0, 0, 0, 0,
                          0, 0, 0, 1, 0, 0,
                          0, 0, 0, 0, 1, 0), nrow=4, ncol=6, byrow = TRUE),
  r_Beta       = matrix(c(1, 0, 1, 0), nrow = 4, ncol = 1),
  plotRoots    = 0 # Don't create plots for tests.
)
m2r2 <- FCVARestn(x, k = 2, r = 2, opt2)


#----------------------------------------------------------------------
# Try variants on restrictions on beta.
#----------------------------------------------------------------------


opt <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  plotRoots    = 0 # Don't create plots for tests.
)
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
m1 <- FCVARestn(x, k = 2, r = 1, opt)



# Or start a new object with otherwise same restrictions.
opt1 <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  R_Beta       = matrix(c(1, 0, 0), nrow = 1, ncol = 3),
  plotRoots    = 0 # Don't create plots for tests.
)
m1r1_0 <- FCVARestn(x, k = 2, r = 1, opt1)



# Or start a new object with otherwise same restrictions.
opt1 <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  R_Beta       = matrix(c(1, 0, 0), nrow = 1, ncol = 3),
  r_Beta       = matrix(c(1), nrow = 1, ncol = 1),
  plotRoots    = 0 # Don't create plots for tests.
)
m1r1_1 <- FCVARestn(x, k = 2, r = 1, opt1)




# Or start a new object with otherwise same restrictions.
opt1 <- FCVARoptions(
  gridSearch   = 0, # Disable grid search in optimization.
  dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
  dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
  constrained  = 0, # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
  R_Beta       = matrix(c(1, 0, 0), nrow = 1, ncol = 3),
  r_Beta       = matrix(c(1.1), nrow = 1, ncol = 1),
  plotRoots    = 0 # Don't create plots for tests.
)
m1r1_11 <- FCVARestn(x, k = 2, r = 1, opt1)



# These models are observationally equivalent,
# with a scalar multiplying and dividing alpha and beta.
# Now try testing hypotheses to throw errors.

# This is a valid test, with beta_1 = 0.
H_m1r1_0 <- FCVARhypoTest(modelUNR = m1, modelR = m1r1_0)

# These are not:
H_m1r1_1 <- FCVARhypoTest(modelUNR = m1, modelR = m1r1_1)
H_m1r1_11 <- FCVARhypoTest(modelUNR = m1, modelR = m1r1_11)



