
# Results from check 7:

# Last result (after some changes):
1 error x | 0 warnings √ | 1 note x

# This time:

> devtools::check()
Updating FCVAR documentation
Loading FCVAR
Writing NAMESPACE
Writing NAMESPACE
Writing plot.FCVAR_grid.Rd
-- Building ----------------------------------------------------------- FCVAR --
  Setting env vars:
  * CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
* CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
* CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
--------------------------------------------------------------------------------
  √  checking for file 'C:\Users\le279259\Documents\Research\FCVAR\GitRepo\FCVAR/DESCRIPTION' ...
-  preparing 'FCVAR': (14.4s)
√  checking DESCRIPTION meta-information ...
-  checking for LF line-endings in source and make files and shell scripts (1.8s)
-  checking for empty or unneeded directories
Removed empty directory 'FCVAR/vignettes'
-  building 'FCVAR_0.0.0.9000.tar.gz'

-- Checking ----------------------------------------------------------- FCVAR --
  Setting env vars:
  * _R_CHECK_CRAN_INCOMING_USE_ASPELL_: TRUE
* _R_CHECK_CRAN_INCOMING_REMOTE_    : FALSE
* _R_CHECK_CRAN_INCOMING_           : FALSE
* _R_CHECK_FORCE_SUGGESTS_          : FALSE
* NOT_CRAN                          : true
-- R CMD check -------------------------------------------------------------------------------------
  -  using log directory 'C:/Users/le279259/AppData/Local/Temp/1/RtmpI7839i/FCVAR.Rcheck'
-  using R version 4.0.2 (2020-06-22)
-  using platform: x86_64-w64-mingw32 (64-bit)
-  using session charset: ISO8859-1
-  using options '--no-manual --as-cran'
√  checking for file 'FCVAR/DESCRIPTION' ...
-  this is package 'FCVAR' version '0.0.0.9000'
-  package encoding: UTF-8
√  checking package namespace information ...
√  checking package dependencies (977ms)
√  checking if this is a source package ...
√  checking if there is a namespace
√  checking for executable files (415ms)
√  checking for hidden files and directories ...
√  checking for portable file names ...
√  checking whether package 'FCVAR' can be installed (4.2s)
√  checking installed package size ...
√  checking package directory
N  checking for future file timestamps (419ms)
unable to verify current time
√  checking DESCRIPTION meta-information (340ms)
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
√  checking S3 generic/method consistency (512ms)
√  checking replacement functions ...
√  checking foreign function calls ...
√  checking R code for possible problems (6.4s)
√  checking Rd files (414ms)
√  checking Rd metadata ...
√  checking Rd line widths ...
√  checking Rd cross-references ...
√  checking for missing documentation entries ...
√  checking for code/documentation mismatches (625ms)
√  checking Rd \usage sections (722ms)
√  checking Rd contents ...
√  checking for unstated dependencies in examples (412ms)
√  checking contents of 'data' directory ...
√  checking data for non-ASCII characters ...
√  checking data for ASCII and uncompressed saves ...
√  checking examples (2h 39m 56.5s)
Examples with CPU (user + system) or elapsed time > 5s
user system elapsed
plot.FCVAR_grid     4714.80   0.41 4719.08
FCVARlikeGrid       4679.27   0.27 4682.30
FCVARestn             29.80   0.00   29.94
FCVARhypoTest         29.80   0.00   29.95
FCVARboot             27.34   0.03   27.38
FCVARbootRank         20.97   0.00   20.98
FCVARforecast          8.12   0.01    8.19
FCVARsimBS             7.56   0.00    7.58
FCVARsim               7.55   0.00    7.58
summary.MVWN_stats     7.36   0.03    7.39
summary.FCVAR_roots    7.31   0.00    7.34
plot.FCVAR_roots       7.30   0.00    7.34
summary.FCVAR_model    7.25   0.01    7.28
MVWNtest               7.25   0.00    7.25
GetCharPolyRoots       7.06   0.01    7.10
summary.FCVAR_ranks    5.36   0.02    5.38
FCVARrankTests         5.22   0.00    5.22
√  checking for unstated dependencies in 'tests' ...
-  checking tests ...
√  Running 'testthat.R' [53s] (53s)
√  checking for non-standard things in the check directory (53s)
√  checking for detritus in the temp directory

See
'C:/Users/le279259/AppData/Local/Temp/1/RtmpI7839i/FCVAR.Rcheck/00check.log'
for details.


-- R CMD check results ------------------------------------------------------- FCVAR 0.0.0.9000 ----
  Duration: 2h 41m 10.8s

> checking for future file timestamps ... NOTE
unable to verify current time

0 errors √ | 0 warnings √ | 1 note x

# One note to correct.
>
