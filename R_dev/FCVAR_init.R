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
# Copied from explorer:
# 'C:\Users\le279259\Documents\Research\FCVAR\GitRepo\FCVAR'
fcvar_dir <- 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR'
setwd(fcvar_dir)

# Since folder already exists, use setup() instead of create().
devtools::setup(path = fcvar_dir)
# Error: 'setup' is not an exported object from 'namespace:devtools'
setup(path = fcvar_dir)

# devtools is no longer current. It was split into several pieces. 
# There is a separate package for package setup. 
library(usethis)
usethis::setup(path = fcvar_dir)
# Does not exist. See below with create_package().

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

# Ok now create the package.

# New project 'FCVAR' is nested inside an existing project 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/', which is rarely a good idea.
# Do you want to create anyway?
#   
#   1: Absolutely not
# 2: No
# 3: Yes
# 
# Selection: Y
# Enter an item from the menu, or 0 to exit
# Selection: Yes
# check Setting active project to 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR'
# check Creating 'R/'
# check Writing 'DESCRIPTION'
# Package: FCVAR
# Title: What the Package Does (One Line, Title Case)
# Version: 0.0.0.9000
# Authors@R (parsed):
#   * First Last <first.last@example.com> [aut, cre] (<https://orcid.org/YOUR-ORCID-ID>)
# Description: What the package does (one paragraph).
# License: What license it uses
# Encoding: UTF-8
# LazyData: true
# check Writing 'NAMESPACE'
# check Writing 'FCVAR.Rproj'
# check Adding '.Rproj.user' to '.gitignore'
# check Adding '^FCVAR\\.Rproj$', '^\\.Rproj\\.user$' to '.Rbuildignore'
# check Opening 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR/' in new RStudio session
# check Setting active project to '<no active project>'



##################################################
# Create documantation

##################################################

# Next steps:
# Edit field in the description
# Create a function in the R folder.
# Add roxygen comments before the function. 

# Use roxygen to build documentation. 
devtools::document()
# Rinse and repeat. 


# Set some folders to be ignored by R build.
devtools::use_build_ignore(c("stata", "R_dev", "MATLAB"))
# Except that is deprecated too.
usethis::use_build_ignore(c("stata", "R_dev", "MATLAB"))
# check Setting active project to 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR'
# check Adding '^stata$', '^R_dev$', '^MATLAB$' to '.Rbuildignore'



# Generate a pdf manual once the documentatino is complete.
devtools::build_manual()
# Hmm ... looks like a package
# Creating pdf output from LaTeX ...
# Saving output to 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR_0.0.0.9000.pdf' ...
# Done


# Generate a sketch of a vignette.
usethis::use_vignette("FCVAR")
# check Setting active project to 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR'
# check Adding 'rmarkdown' to Suggests field in DESCRIPTION
# check Writing 'vignettes/FCVAR.Rmd'
# dot Modify 'vignettes/FCVAR.Rmd'

# Generate a sketch of an article.
usethis::use_article("FCVAR")
# Overwrite pre-existing file 'vignettes/FCVAR.Rmd'?
#   
#   1: Yeah
# 2: Negative
# 3: Not now
# 
# Selection: 1
# check Writing 'vignettes/FCVAR.Rmd'
# dot Modify 'vignettes/FCVAR.Rmd'
# check Adding '^vignettes/FCVAR\\.Rmd$' to '.Rbuildignore'

# Underwhelming. Only creates a folder. 

# Try to build:
devtools::build_vignettes()
# NULL
# No actions that I can detect with git status.



# Iterations on R package.
# load_all() is the key step in this "lather, rinse, repeat" 
# cycle of package development:
#   
# 1.  Tweak a function definition.
# 2. 
load_all()
# 3. Try out the change by running a small example or some tests.


# Examples

# Data for examples. 
votingJNP2014 <- read.csv('data_JNP2014.csv')
colnames(votingJNP2014)
# [1] "lib"    "pc"     "ir_can" "ir_us"  "un_can" "un_us"
# x1 is the only one used.
colnames(votingJNP2014)[c(1, 3, 5)]
x1 <- votingJNP2014[, c("lib", "ir_can", "un_can")]

# Then run
# usethis::use_data() 
usethis::use_data(votingJNP2014, votingJNP2014)
# Warning: Saving duplicates only once: 'votingJNP2014'
# check Creating 'data/'
# check Saving 'votingJNP2014' to 'data/votingJNP2014.rda'

# You should also make sure that the data has been optimally compressed:
# Run 
tools::checkRdaFiles() 
# to determine the best compression for each file.

# Re-run 
usethis::use_data() 
# with compress set to that optimal value. 
# If you've lost the code for recreating the files, you can use 
tools::resaveRdaFiles() 
# to re-save in place.




