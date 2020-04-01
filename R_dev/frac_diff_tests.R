################################################################################
# 
# Testing algorithms for fast fractional differencing
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
# Testing algorithms for fast fractional differencing:
#   1. Simple time domain calculation with a filter
#   2. Fast fractional differencing from Jensen and Nielsen
#   3. Fast fractional differencing from the fracdiff R package
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
# source('EstOptions.R')
# Contains a single function EstOptions() with default options. 

# Load library with Main Estimation Function for FCVAR
# source('FCVAR_estn.R')

# Load library of Functions for FCVAR Estimation
source('FCVAR_lower.R')

# Load library of Pre- and Post-Estimation Functions for FCVAR
# source('FCVAR_higher.R')

# Load function for fast fractional differencing by Jensen and Nielsen.
source('fracdiff.R')


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
# Comparison with simple series.
################################################################################

d <- 0.2
T_simple <- 100
# test_series <- matrix(rep(1, T_simple), nrow = T_simple, ncol = 1) # FD is zero
test_series <- matrix(rnorm(T_simple), nrow = T_simple, ncol = 1)
# test_series <- matrix(data[, 1])
# test_series <- matrix(data[, 2])
# test_series <- matrix(data[, 3]) # Big difference (at beginning of series, bigger with small d)
# test_series <- matrix(data[, 5]) # Big difference (at beginning of series)


comp_simple <- data.frame(matrix(NA, nrow = length(test_series), ncol = 3))
colnames(comp_simple) <- c('time_domain', 'JN_FFT', 'FD_FFT')


source('FCVAR_lower.R')
source('fracdiff.R')


comp_simple[, 'time_domain'] <- FracDiff_filter(test_series, d)
comp_simple[, 'JN_FFT'] <- fracdiff_JN(test_series, d)
comp_simple[, 'FD_FFT'] <- diffseries(test_series, d)


summary(comp_simple)
head(comp_simple)
tail(comp_simple)

plot(comp_simple[, 'JN_FFT'], comp_simple[, 'time_domain'])
plot(comp_simple[, 'JN_FFT'], comp_simple[, 'FD_FFT'])

sum(abs(comp_simple[, 'JN_FFT'] - comp_simple[, 'time_domain']))
sum(abs(comp_simple[, 'JN_FFT'] - comp_simple[, 'FD_FFT']))


sum(abs(comp_simple[, 'JN_FFT'] - comp_simple[, 'FD_FFT']) > 10^(-1))
which(abs(comp_simple[, 'JN_FFT'] - comp_simple[, 'FD_FFT']) > 10^(-1))



################################################################################
# Test the multivariate versions.
################################################################################

# d <- 0.2
d <- 0.8

source('FCVAR_lower.R')
source('fracdiff.R')

test_mat <- data[, c(3, 5)]

comp_multiple <- data.frame(matrix(NA, nrow = nrow(test_mat), ncol = 6))
colnames(comp_multiple) <- c('time_domain_1', 'FCVAR_FFT_1', 'JN_FFT_1', 
                                     'time_domain_2', 'FCVAR_FFT_2', 'JN_FFT_2')


comp_multiple[, c('time_domain_1', 'time_domain_2')] <- FracDiff_filter(test_mat, d)
comp_multiple[, c('FCVAR_FFT_1', 'FCVAR_FFT_2')] <- FracDiff(test_mat, d)
comp_multiple[, c('JN_FFT_1', 'JN_FFT_2')] <- fracdiff_JN(test_mat, d)


# Compare series.
summary(comp_multiple[, c(1,2,3)])
summary(comp_multiple[, c(4,5,6)])



sum(abs(comp_multiple[, 'time_domain_1'] - comp_multiple[, 'FCVAR_FFT_1']))
sum(abs(comp_multiple[, 'time_domain_2'] - comp_multiple[, 'FCVAR_FFT_2']))


  

################################################################################
# End
################################################################################
