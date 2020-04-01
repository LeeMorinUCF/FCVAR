################################################################################
# 
# Creating the FCVAR Package
# 
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 29, 2020
# 
################################################################################
# 
# This program collects the tools needed to construct
#   an R package for the FCVAR model.
#   It uses the packages in Hadley Wickham's book "R Packages"
# 
################################################################################


################################################################################
# Clearing Workspace and Declaring Packages
################################################################################

# Clear workspace.
rm(list = ls(all = TRUE))

# Set working directory.
wd_path <- '~/Research/FCVAR/GitRepo/FCVAR/R'
setwd(wd_path)

# Load package for development tools.
# install.packages('devtools')
library(devtools)

# Load package for generating R documentation.
# install.packages('roxygen2')
library(roxygen2)

# Load package for testing functions.
# install.packages('testthat')
library(testthat)


# Load package for generating R documentation
# with snippets of R code.
# install.packages('knitr')
library(knitr)


# Check that R version is sufficient.
library(rstudioapi)
rstudioapi::isAvailable("0.99.149")
# [1] TRUE


# Check that everything is installed:
has_devel()
# Need grown-up permission to install. 
# Diagnosing the problem.
pkgbuild::check_build_tools(debug = TRUE)
# Have you ever seen the movie 'Inception'?



devtools::session_info()


# Also test Rcpp is installed correctly. 
library(Rcpp)

# Rcpp example.
sourceCpp('Rcpp_tests/cumsum1_test.cpp')

print('Testing with c(1, 2, 3, 4, 5):')
# Rccp version:
cumsum1(c(1.0, 2.0, 3.0, 4.0, 5.0))
# Baes R function:
cumsum(c(1.0, 2.0, 3.0, 4.0, 5.0))


# Simpler example:
sourceCpp('Rcpp_tests/timesTwo.cpp')

# Test
timesTwo(21)



################################################################################
# End
################################################################################

