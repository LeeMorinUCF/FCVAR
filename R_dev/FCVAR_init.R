##################################################
# 
# aggregress Scratchpad
# 
# Lealand Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# February 15, 2019
# 
##################################################
# 
# A set of commands to initialize the FCVAR package
# 
# Dependencies:
#   None.
# 
##################################################


##################################################
# Preparing the Workspace
##################################################

# Clear workspace.
rm(list=ls(all=TRUE))




# Check that the package development tools are installed correctly
library(devtools)
# Loading required package: usethis
# Warning messages:
#   1: package 'devtools' was built under R version 3.5.3 
#   2: package 'usethis' was built under R version 3.5.3 
has_devel()
# Your system is ready to build packages!



##################################################
# Package initialization. 
##################################################

fcvar_dir <- '~/Research/FCVAR/GitRepo/FCVAR'
setwd(fcvar_dir)

# Since folder already exists, use setup() instead of create().
devtools::setup(path = fcvar_dir)
# Error: 'setup' is not an exported object from 'namespace:devtools'
setup(path = fcvar_dir)

# devtools is no longer current. It was split into several pieces. 
# There is a separate package for package setup. 
library(usethis)
usethis::setup(path = fcvar_dir)

# Check that the name is available.
install.packages('available')
library(available)
available('FCVAR')
# Urban Dictionary can contain potentially offensive results,
# should they be included? [Y]es / [N]o:
#   1: Y
# -- FCVAR -----------------------------------------------------------------------
#   Name valid: check
# Available on CRAN: check
# Available on Bioconductor: check
# Available on GitHub:  check
# Abbreviations: http://www.abbreviations.com/FCVAR
# Wikipedia: https://en.wikipedia.org/wiki/FCVAR
# Wiktionary: https://en.wiktionary.org/wiki/FCVAR
# Urban Dictionary:
#   Not found.
# Sentiment:???
#   Warning message:
#   package 'tidytext' was built under R version 3.5.3 

# This means we are G2G. 



