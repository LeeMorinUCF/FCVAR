##################################################
#
# FCVAR Initialization Scratchpad
#
# Lealand Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
#
# March 30, 2020
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
usethis::use_build_ignore(c("stata", "R_dev", "MATLAB", "article"))



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
tools::checkRdaFiles('data/votingJNP2014.rda')
# size ASCII compress version
# data/votingJNP2014.rda 9233 FALSE    bzip2       2
# to determine the best compression for each file.
# So, the best compression is 'bzip2', which is the default, so we're good.
# The file is already quite small.

# Otherwise:
# Re-run
usethis::use_data(votingJNP2014, votingJNP2014, compress = 'whatever')
# with compress set to that optimal value.
# If you've lost the code for recreating the files, you can use
tools::resaveRdaFiles()
# to re-save in place.



# Testing
usethis::use_testthat()
# check Adding 'testthat' to Suggests field in DESCRIPTION
# check Creating 'tests/testthat/'
# check Writing 'tests/testthat.R'
# dot Call `use_test()` to initialize a basic test file and open it for editing.


# Dependencies.
# Depends on the fracdist package on GitHUb.
# Instructions from
# https://stackoverflow.com/questions/30493388/create-an-r-package-that-depends-on-another-r-package-located-on-github
#
#
# Package developer POV:
#
#   1) do:
#
#   devtools::use_package("jvamisc")
#   devtools::document()
# to add the dependency in the Imports field of your DESCRIPTION file.
#
# 2) manually add a field "Remotes:" in the DESCRIPTION file, specifying where on github R should look for the package:
#
#   #in DESCRIPTION
#   Imports: ...,
#     jvamisc,
#     ...
# Remotes: JVAdams/jvamisc

# So, here we go:
# devtools::use_package("fracdist")
# Error: 'use_package' is not an exported object from 'namespace:devtools'
# usethis::use_package('fracdist')
# check Setting active project to 'C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR'
# check Adding 'fracdist' to Imports field in DESCRIPTION
# - Refer to functions with `fracdist::fun()`
# Added to the bottom of DESCRIPTION:
# Remotes: LeeMorinUCF/fracdist





# Checking status of package:

# devtools::check()
# Updating FCVAR documentation
# Writing NAMESPACE
# Loading FCVAR
# Writing NAMESPACE
# -- Building ----------------------------------------------------------- FCVAR --
#   Setting env vars:
#   * CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
# * CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
# * CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
# --------------------------------------------------------------------------------
#   check  checking for file 'C:\Users\le279259\Documents\Research\FCVAR\GitRepo\FCVAR/DESCRIPTION' (605ms)
# -  preparing 'FCVAR': (18.1s)
# check  checking DESCRIPTION meta-information ...
# -  checking for LF line-endings in source and make files and shell scripts (1.8s)
# -  checking for empty or unneeded directories
# Removed empty directory 'FCVAR/vignettes'
# -  looking to see if a 'data/datalist' file should be added
# -  building 'FCVAR_0.0.0.9000.tar.gz'
#
# -- Checking ----------------------------------------------------------- FCVAR --
#   Setting env vars:
#   * _R_CHECK_CRAN_INCOMING_USE_ASPELL_: TRUE
# * _R_CHECK_CRAN_INCOMING_REMOTE_    : FALSE
# * _R_CHECK_CRAN_INCOMING_           : FALSE
# * _R_CHECK_FORCE_SUGGESTS_          : FALSE
# * NOT_CRAN                          : true
# -- R CMD check ---------------------------------------------------------------------
#   -  using log directory 'C:/Users/le279259/AppData/Local/Temp/1/RtmpuMbWJI/FCVAR.Rcheck'
# -  using R version 3.5.1 (2018-07-02)
# -  using platform: x86_64-w64-mingw32 (64-bit)
# -  using session charset: ISO8859-1
# -  using options '--no-manual --as-cran' (470ms)
# check  checking for file 'FCVAR/DESCRIPTION'
# -  this is package 'FCVAR' version '0.0.0.9000'
# -  package encoding: UTF-8
# check  checking package namespace information ...
# check  checking package dependencies (1.1s)
# check  checking if this is a source package ...
# check  checking if there is a namespace
# check  checking for executable files (642ms)
# check  checking for hidden files and directories ...
# check  checking for portable file names ...
# check  checking serialization versions ...
# W  checking whether package 'FCVAR' can be installed (4.8s)
# Found the following significant warnings:
#   Note: possible error in 'matrix(1, nrow = T, nol = 1)': unused argument (nol = 1)
# See 'C:/Users/le279259/AppData/Local/Temp/1/RtmpuMbWJI/FCVAR.Rcheck/00install.out' for details.
# Information on the location(s) of code generating the 'Note's can be
# obtained by re-running with environment variable R_KEEP_PKG_SOURCE set
# to 'yes'.
# check  checking installed package size ...
# check  checking package directory
# check  checking DESCRIPTION meta-information (350ms)
# N  checking top-level files
# Non-standard file/directory found at top level:
#   'article'
# check  checking for left-over files ...
# W  checking index information ...
# Demos with missing or empty index information:
#   FCVAR_extension_JNP2014
# FCVAR_replication_JNP2014
# See sections 'The INDEX file' and 'Package subdirectories' in the
# 'Writing R Extensions' manual.
# check  checking package subdirectories ...
# check  checking R files for non-ASCII characters ...
# check  checking R files for syntax errors ...
# check  checking whether the package can be loaded ...
# check  checking whether the package can be loaded with stated dependencies ...
# check  checking whether the package can be unloaded cleanly ...
# check  checking whether the namespace can be loaded with stated dependencies ...
# check  checking whether the namespace can be unloaded cleanly ...
# check  checking loading without being on the library search path ...
# check  checking dependencies in R code ...
# W  checking S3 generic/method consistency (514ms)
# print:
#   function(x, ...)
#     print.FCVARestn:
#   function(results, k, r, p, T, opt)
#
#     print:
#   function(x, ...)
#     print.LagSelect:
#   function(stats, kmax, r, p, T, order, opt)
#
#     print:
#   function(x, ...)
#     print.RankTests:
#   function(stats, k, p, T, opt)
#
#     See section 'Generic functions and methods' in the 'Writing R
#    Extensions' manual.
# check  checking replacement functions ...
# check  checking foreign function calls ...
# N  checking R code for possible problems (6.1s)
# FCVARestn: no visible global function definition for 'optim'
# FCVARestn: no visible global function definition for 'constrOptim'
# FCVARoptionUpdates: no visible binding for global variable 'end'
# FCVARsim: no visible global function definition for 'rnorm'
# FCVARsimBS: no visible global function definition for 'rnorm'
# FCVARsimBS: possible error in matrix(1, nrow = T, nol = 1): unused
# argument (nol = 1)
# FracDiff: no visible global function definition for 'nextn'
# FracDiff: no visible global function definition for 'fft'
# HypoTest: no visible global function definition for 'pchisq'
# LMtest: no visible global function definition for 'pchisq'
# LagSelect: no visible global function definition for 'pchisq'
# LikeGridSearch: no visible global function definition for 'optim'
# Qtest: no visible global function definition for 'pchisq'
# plot.GetCharPolyRoots: no visible global function definition for 'pdf'
# plot.GetCharPolyRoots: no visible global function definition for 'png'
# plot.GetCharPolyRoots: no visible global function definition for 'plot'
# plot.GetCharPolyRoots: no visible global function definition for
# 'lines'
# plot.GetCharPolyRoots: no visible global function definition for
# 'points'
# plot.GetCharPolyRoots: no visible global function definition for
# 'dev.off'
# plot.LikeGridSearch: no visible global function definition for 'pdf'
# plot.LikeGridSearch: no visible global function definition for 'png'
# plot.LikeGridSearch: no visible global function definition for
# 'rainbow'
# plot.LikeGridSearch: no visible global function definition for 'par'
# plot.LikeGridSearch: no visible global function definition for 'persp'
# plot.LikeGridSearch: no visible global function definition for 'plot'
# plot.LikeGridSearch: no visible global function definition for
# 'dev.off'
# print.MVWNtest: no visible binding for global variable 'p'
# Undefined global functions or variables:
#   constrOptim dev.off end fft lines nextn optim p par pchisq pdf persp
# plot png points rainbow rnorm
# Consider adding
# importFrom("grDevices", "dev.off", "pdf", "png", "rainbow")
# importFrom("graphics", "lines", "par", "persp", "plot", "points")
# importFrom("stats", "constrOptim", "end", "fft", "nextn", "optim",
#            "pchisq", "rnorm")
# to your NAMESPACE file.
# check  checking Rd files (507ms)
# check  checking Rd metadata ...
# check  checking Rd line widths ...
# check  checking Rd cross-references (413ms)
# check  checking for missing documentation entries ...
# check  checking for code/documentation mismatches (734ms)
# W  checking Rd \usage sections (414ms)
# Undocumented arguments in documentation object 'FCVARlikeFull'
# 'beta' 'rho'
# Documented arguments not in \usage in documentation object 'FCVARlikeFull':
#   'betaHat' 'rhoHat'
#
# Undocumented arguments in documentation object 'GetFreeParams'
# 'p'
# Documented arguments not in \usage in documentation object 'GetFreeParams':
#   'x'
#
# Functions with \usage entries need to have the appropriate \alias
# entries, and all their arguments documented.
# The \usage entries must correspond to syntactically valid R code.
# See chapter 'Writing R documentation files' in the 'Writing R
#    Extensions' manual.
# check  checking Rd contents (721ms)
# check  checking for unstated dependencies in examples (613ms)
# check  checking contents of 'data' directory
# check  checking data for non-ASCII characters ...
# W  checking data for ASCII and uncompressed saves


# Checking examples...

# Fix the above and cut the examples short.
# i.e. B = 9 vs 999 in comments.

