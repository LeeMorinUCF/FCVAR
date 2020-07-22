

# Results from check 5:

# Last result (the never-ending errors):
1 error x | 1 warning x | 0 notes √

# This time:

> devtools::check()
Updating FCVAR documentation
Writing NAMESPACE
Loading FCVAR
Writing NAMESPACE
Writing print.MVWNtest.Rd
-- Building ----------------------------------------------------------- FCVAR --
  Setting env vars:
  * CFLAGS    : -Wall -pedantic
* CXXFLAGS  : -Wall -pedantic
* CXX11FLAGS: -Wall -pedantic
--------------------------------------------------------------------------------
  √  checking for file 'C:\Users\le279259\Documents\Research\FCVAR\GitRepo\FCVAR/DESCRIPTION' ...
-  preparing 'FCVAR': (13.3s)
√  checking DESCRIPTION meta-information ...
-  checking for LF line-endings in source and make files and shell scripts (1.7s)
-  checking for empty or unneeded directories
Removed empty directory 'FCVAR/vignettes'
-  looking to see if a 'data/datalist' file should be added
-  building 'FCVAR_0.0.0.9000.tar.gz'

-- Checking ----------------------------------------------------------- FCVAR --
  Setting env vars:
  * _R_CHECK_CRAN_INCOMING_USE_ASPELL_: TRUE
* _R_CHECK_CRAN_INCOMING_REMOTE_    : FALSE
* _R_CHECK_CRAN_INCOMING_           : FALSE
* _R_CHECK_FORCE_SUGGESTS_          : FALSE
* NOT_CRAN                          : true
-- R CMD check -------------------------------------------------------------------
  -  using log directory 'C:/Users/le279259/AppData/Local/Temp/1/RtmpeWh3xu/FCVAR.Rcheck'
-  using R version 3.5.1 (2018-07-02)
-  using platform: x86_64-w64-mingw32 (64-bit)
-  using session charset: ISO8859-1
-  using options '--no-manual --as-cran' (359ms)
√  checking for file 'FCVAR/DESCRIPTION'
-  this is package 'FCVAR' version '0.0.0.9000'
-  package encoding: UTF-8
√  checking package namespace information ...
√  checking package dependencies (905ms)
√  checking if this is a source package ...
√  checking if there is a namespace
√  checking for executable files (416ms)
√  checking for hidden files and directories ...
√  checking for portable file names ...
√  checking whether package 'FCVAR' can be installed (4.4s)
√  checking installed package size ...
√  checking package directory
√  checking DESCRIPTION meta-information (339ms)
√  checking top-level files
√  checking for left-over files ...
√  checking index information ...
√  checking package subdirectories ...
√  checking R files for non-ASCII characters ...
√  checking R files for syntax errors ...
√  checking whether the package can be loaded ...
√  checking whether the package can be loaded with stated dependencies ...
√  checking whether the package can be unloaded cleanly ...
√  checking whether the namespace can be loaded with stated dependencies ...
√  checking whether the namespace can be unloaded cleanly ...
√  checking loading without being on the library search path ...
√  checking dependencies in R code ...
W  checking S3 generic/method consistency (515ms)
print:
  function(x, ...)
    print.FCVARestn:
  function(results, k, r, p, T, opt)

    print:
  function(x, ...)
    print.FCVARlagSelect:
  function(stats, kmax, r, p, T, order, opt)

    print:
  function(x, ...)
    print.FCVARrankTests:
  function(stats, k, p, T, opt)

    print:
  function(x, ...)
    print.GetCharPolyRoots:
  function(cPolyRoots)

    print:
  function(x, ...)
    print.MVWNtest:
  function(stats, maxlag, p)

    plot:
  function(x, ...)
    plot.FCVARlikeGrid:
  function(likeGrid_params, k, r, opt, file, file_ext, main)

    plot:
  function(x, ...)
    plot.GetCharPolyRoots:
  function(cPolyRoots, b, file, file_ext, xlim, ylim, main)

    See section 'Generic functions and methods' in the 'Writing R
   Extensions' manual.
√  checking replacement functions ...
√  checking foreign function calls ...
√  checking R code for possible problems (5.7s)
√  checking Rd files (409ms)
√  checking Rd metadata ...
√  checking Rd line widths ...
√  checking Rd cross-references ...
√  checking for missing documentation entries ...
√  checking for code/documentation mismatches (731ms)
√  checking Rd \usage sections (824ms)
√  checking Rd contents ...
√  checking for unstated dependencies in examples (615ms)
√  checking contents of 'data' directory
√  checking data for non-ASCII characters ...
√  checking data for ASCII and uncompressed saves ...
√  checking examples (7m 7.8s)
Examples with CPU or elapsed time > 5s
user system elapsed
plot.FCVARlikeGrid     127.48   0.01  127.59
FCVARhypoTest           28.24   0.00   28.33
FCVARestn               28.14   0.00   28.29
FCVARboot               26.50   0.03   26.53
FCVARlikeGrid           24.89   0.02   24.93
FCVARbootRank           20.17   0.00   20.17
FCVARhess               12.54   0.02   12.59
FCVARsimBS               8.11   0.03    8.15
FCVARsim                 7.92   0.03    8.00
FCVARforecast            7.79   0.00    7.81
SEmat2vecU               7.59   0.02    7.63
LMtest                   7.53   0.00    7.56
GetCharPolyRoots         7.48   0.00    7.56
print.MVWNtest           7.44   0.00    7.50
Qtest                    7.41   0.00    7.47
MVWNtest                 7.35   0.00    7.38
GetResiduals             7.34   0.00    7.39
Lbk                      7.16   0.01    7.22
print.GetCharPolyRoots   7.17   0.00    7.26
TransformData            7.14   0.00    7.17
GetEstimates             7.11   0.02    7.16
FCVARlike                7.03   0.00    7.04
FCVARlikeFull            7.03   0.00    7.04
SEvec2matU               6.97   0.00    7.00
print.FCVARestn          6.81   0.02    6.86
plot.GetCharPolyRoots    6.77   0.00    6.81
FCVARrankTests           5.54   0.00    5.58
print.FCVARrankTests     5.12   0.00    5.12
√  checking for unstated dependencies in 'tests' ...
-  checking tests ...
√  Running 'testthat.R' [54s] (54.4s)

See
'C:/Users/le279259/AppData/Local/Temp/1/RtmpeWh3xu/FCVAR.Rcheck/00check.log'
for details.


-- R CMD check results ------------------------------------- FCVAR 0.0.0.9000 ----
  Duration: 8m 23s

> checking S3 generic/method consistency ... WARNING
print:
  function(x, ...)
    print.FCVARestn:
  function(results, k, r, p, T, opt)

    print:
  function(x, ...)
    print.FCVARlagSelect:
  function(stats, kmax, r, p, T, order, opt)

    print:
  function(x, ...)
    print.FCVARrankTests:
  function(stats, k, p, T, opt)

    print:
  function(x, ...)
    print.GetCharPolyRoots:
  function(cPolyRoots)

    print:
  function(x, ...)
    print.MVWNtest:
  function(stats, maxlag, p)

    plot:
  function(x, ...)
    plot.FCVARlikeGrid:
  function(likeGrid_params, k, r, opt, file, file_ext, main)

    plot:
  function(x, ...)
    plot.GetCharPolyRoots:
  function(cPolyRoots, b, file, file_ext, xlim, ylim, main)

    See section 'Generic functions and methods' in the 'Writing R
  Extensions' manual.

0 errors √ | 1 warning x | 0 notes √
>
