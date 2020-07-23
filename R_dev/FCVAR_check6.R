
# Results from check 6:

# Last result:
0 errors √ | 1 warning x | 0 notes √

# This time:

> devtools::check()
Updating FCVAR documentation
Writing NAMESPACE
Loading FCVAR
Writing NAMESPACE
Writing plot.FCVAR_roots.Rd
-- Building ----------------------------------------------------------- FCVAR --
  Setting env vars:
  * CFLAGS    : -Wall -pedantic
* CXXFLAGS  : -Wall -pedantic
* CXX11FLAGS: -Wall -pedantic
--------------------------------------------------------------------------------
  √  checking for file 'C:\Users\le279259\Documents\Research\FCVAR\GitRepo\FCVAR/DESCRIPTION' ...
-  preparing 'FCVAR': (12.6s)
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
-- R CMD check ----------------------------------------------------------------------------------------------
  -  using log directory 'C:/Users/le279259/AppData/Local/Temp/1/RtmpAhLgcJ/FCVAR.Rcheck'
-  using R version 3.5.1 (2018-07-02)
-  using platform: x86_64-w64-mingw32 (64-bit)
-  using session charset: ISO8859-1
-  using options '--no-manual --as-cran'
√  checking for file 'FCVAR/DESCRIPTION' ...
-  this is package 'FCVAR' version '0.0.0.9000'
-  package encoding: UTF-8
√  checking package namespace information ...
√  checking package dependencies (929ms)
√  checking if this is a source package ...
√  checking if there is a namespace
√  checking for executable files (466ms)
√  checking for hidden files and directories ...
√  checking for portable file names ...
√  checking whether package 'FCVAR' can be installed (3.7s)
√  checking installed package size ...
√  checking package directory
√  checking DESCRIPTION meta-information (359ms)
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
√  checking S3 generic/method consistency (440ms)
√  checking replacement functions ...
√  checking foreign function calls (339ms)
√  checking R code for possible problems (5.2s)
√  checking Rd files (434ms)
√  checking Rd metadata ...
√  checking Rd line widths ...
√  checking Rd cross-references ...
√  checking for missing documentation entries ...
√  checking for code/documentation mismatches (562ms)
√  checking Rd \usage sections (752ms)
√  checking Rd contents ...
√  checking for unstated dependencies in examples (539ms)
√  checking contents of 'data' directory
√  checking data for non-ASCII characters ...
√  checking data for ASCII and uncompressed saves ...
√  checking examples (6m 14.1s)
Examples with CPU or elapsed time > 5s
user system elapsed
plot.FCVAR_grid     114.16   0.01  114.21
FCVARestn            25.06   0.00   25.19
FCVARhypoTest        25.04   0.00   25.17
FCVARboot            23.57   0.00   23.56
FCVARlikeGrid        20.87   0.00   20.89
FCVARbootRank        17.89   0.00   17.90
FCVARhess            11.11   0.00   11.13
FCVARforecast         6.90   0.00    6.96
FCVARlikeFull         6.64   0.01    6.65
FCVARsim              6.55   0.00    6.56
FCVARsimBS            6.50   0.00    6.53
summary.MVWN_stats    6.31   0.00    6.31
Qtest                 6.30   0.00    6.31
LMtest                6.23   0.00    6.25
MVWNtest              6.22   0.00    6.23
plot.FCVAR_roots      6.16   0.02    6.17
summary.FCVAR_roots   6.18   0.00    6.19
Lbk                   6.13   0.03    6.18
GetCharPolyRoots      6.14   0.00    6.16
GetResiduals          6.12   0.00    6.16
SEmat2vecU            6.12   0.00    6.13
TransformData         6.11   0.00    6.11
GetEstimates          6.09   0.02    6.13
summary.FCVAR_model   6.09   0.00    6.14
SEvec2matU            6.08   0.00    6.10
FCVARlike             6.07   0.00    6.08
√  checking for unstated dependencies in 'tests' ...
-  checking tests ...
√  Running 'testthat.R' [45s] (44.9s)


-- R CMD check results ---------------------------------------------------------------- FCVAR 0.0.0.9000 ----
  Duration: 7m 18.1s

0 errors √ | 0 warnings √ | 0 notes √

Mission accomplished!
>
