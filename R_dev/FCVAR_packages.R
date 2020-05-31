##################################################
#
# FCVAR Package Tests
#
# Lealand Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
#
# May 25, 2020
#
##################################################
#
# Testing competing R packages related to the FCVAR
#
# Dependencies:
#   None.
#
##################################################


##################################################
# Packages relating to Cointegration
##################################################

#------------------------------------------------
# The aTSA package
#------------------------------------------------
library(aTSA)

coint.test()


#------------------------------------------------
# The egcm package
#------------------------------------------------
library(egcm)

#------------------------------------------------
# The tseries package
#------------------------------------------------
library(tseries)

#------------------------------------------------
# The cointReg package
#------------------------------------------------
library(cointReg)

#------------------------------------------------
# The urca package
#------------------------------------------------

library(urca)
ca.po()

cajo.test()

#------------------------------------------------
# The tsDyn package
#------------------------------------------------




##################################################
# Packages relating to Long Memory
##################################################

#------------------------------------------------
# The fracdiff package
#------------------------------------------------

library(fracdiff)

x_test <- matrix(rnorm(1000))
x_test <- matrix(cumsum(rnorm(1000)))
d_test <- 0.5

ds_test <- fracdiff::diffseries(x_test, d = d_test)
fd_test <- FCVAR::FracDiff(x_test - mean(x_test), d = d_test)

summary(ds_test - fd_test)

sum(abs(ds_test - fd_test) > 10^(-13))



sum(fracdiff::diffseries(x_test, d = d_test) ==
  FCVAR::FracDiff(x_test - mean(x_test), d = d_test))





#------------------------------------------------
# The arfima package
#------------------------------------------------
library(arfima)

#------------------------------------------------
# The nsarfima package
#------------------------------------------------
library(nsarfima)

#------------------------------------------------
# The LongMemoryTS package
#------------------------------------------------
library(LongMemoryTS)





##################################################
# End
##################################################
