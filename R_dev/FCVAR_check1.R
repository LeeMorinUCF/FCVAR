
# Results from check 1:


> devtools::check()
Updating FCVAR documentation
Writing NAMESPACE
Loading FCVAR
Writing NAMESPACE
-- Building ----------------------------------------------------------- FCVAR --
  Setting env vars:
  * CFLAGS    : -Wall -pedantic
* CXXFLAGS  : -Wall -pedantic
* CXX11FLAGS: -Wall -pedantic
--------------------------------------------------------------------------------
  √  checking for file 'C:\Users\le279259\Documents\Research\FCVAR\GitRepo\FCVAR/DESCRIPTION' ...
-  preparing 'FCVAR': (19s)
√  checking DESCRIPTION meta-information ...
-  checking for LF line-endings in source and make files and shell scripts (1.9s)
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
-  using options '--no-manual --as-cran' (504ms)
√  checking for file 'FCVAR/DESCRIPTION' ...
-  this is package 'FCVAR' version '0.0.0.9000'
-  package encoding: UTF-8
√  checking package namespace information
√  checking package dependencies (1s)
√  checking if this is a source package ...
√  checking if there is a namespace
√  checking for executable files (638ms)
√  checking for hidden files and directories ...
√  checking for portable file names ...
√  checking serialization versions ...
W  checking whether package 'FCVAR' can be installed (5.2s)
Found the following significant warnings:
  Note: possible error in 'matrix(1, nrow = T, nol = 1)': unused argument (nol = 1)
See 'C:/Users/le279259/AppData/Local/Temp/1/RtmpeWh3xu/FCVAR.Rcheck/00install.out' for details.
Information on the location(s) of code generating the 'Note's can be
obtained by re-running with environment variable R_KEEP_PKG_SOURCE set
to 'yes'.
√  checking installed package size ...
√  checking package directory
√  checking DESCRIPTION meta-information (339ms)
N  checking top-level files
Non-standard file/directory found at top level:
  'article'
√  checking for left-over files ...
W  checking index information ...
Demos with missing or empty index information:
  FCVAR_extension_JNP2014
FCVAR_replication_JNP2014
See sections 'The INDEX file' and 'Package subdirectories' in the
'Writing R Extensions' manual.
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
W  checking S3 generic/method consistency (514ms)
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

    See section 'Generic functions and methods' in the 'Writing R
   Extensions' manual.
√  checking replacement functions ...
√  checking foreign function calls ...
N  checking R code for possible problems (6.1s)
FCVARestn: no visible global function definition for 'optim'
FCVARestn: no visible global function definition for 'constrOptim'
FCVARhypoTest: no visible global function definition for 'pchisq'
FCVARlagSelect: no visible global function definition for 'pchisq'
FCVARlikeGrid: no visible global function definition for 'optim'
FCVARoptionUpdates: no visible binding for global variable 'end'
FCVARsim: no visible global function definition for 'rnorm'
FCVARsimBS: no visible global function definition for 'rnorm'
FCVARsimBS: possible error in matrix(1, nrow = T, nol = 1): unused
argument (nol = 1)
FracDiff: no visible global function definition for 'nextn'
FracDiff: no visible global function definition for 'fft'
LMtest: no visible global function definition for 'pchisq'
Qtest: no visible global function definition for 'pchisq'
plot.FCVARlikeGrid: no visible global function definition for 'pdf'
plot.FCVARlikeGrid: no visible global function definition for 'png'
plot.FCVARlikeGrid: no visible global function definition for 'rainbow'
plot.FCVARlikeGrid: no visible global function definition for 'par'
plot.FCVARlikeGrid: no visible global function definition for 'persp'
plot.FCVARlikeGrid: no visible global function definition for 'plot'
plot.FCVARlikeGrid: no visible global function definition for 'dev.off'
plot.GetCharPolyRoots: no visible global function definition for 'pdf'
plot.GetCharPolyRoots: no visible global function definition for 'png'
plot.GetCharPolyRoots: no visible global function definition for 'plot'
plot.GetCharPolyRoots: no visible global function definition for
'lines'
plot.GetCharPolyRoots: no visible global function definition for
'points'
plot.GetCharPolyRoots: no visible global function definition for
'dev.off'
print.MVWNtest: no visible binding for global variable 'p'
Undefined global functions or variables:
  constrOptim dev.off end fft lines nextn optim p par pchisq pdf persp
plot png points rainbow rnorm
Consider adding
importFrom("grDevices", "dev.off", "pdf", "png", "rainbow")
importFrom("graphics", "lines", "par", "persp", "plot", "points")
importFrom("stats", "constrOptim", "end", "fft", "nextn", "optim",
           "pchisq", "rnorm")
to your NAMESPACE file.
√  checking Rd files (507ms)
√  checking Rd metadata ...
√  checking Rd line widths ...
√  checking Rd cross-references (413ms)
√  checking for missing documentation entries ...
√  checking for code/documentation mismatches (734ms)
W  checking Rd \usage sections (414ms)
Undocumented arguments in documentation object 'FCVARlikeFull'
'beta' 'rho'
Documented arguments not in \usage in documentation object 'FCVARlikeFull':
  'betaHat' 'rhoHat'

Undocumented arguments in documentation object 'GetFreeParams'
'p'
Documented arguments not in \usage in documentation object 'GetFreeParams':
  'x'

Functions with \usage entries need to have the appropriate \alias
entries, and all their arguments documented.
The \usage entries must correspond to syntactically valid R code.
See chapter 'Writing R documentation files' in the 'Writing R
   Extensions' manual.
√  checking Rd contents (724ms)
√  checking for unstated dependencies in examples (613ms)
√  checking contents of 'data' directory
√  checking data for non-ASCII characters ...
W  checking data for ASCII and uncompressed saves ...
Warning: package needs dependence on R (>= 2.10)
E  checking examples (24.8s)
Running examples in 'FCVAR-Ex.R' failed
The error most likely occurred in:

  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: FCVARbootRank
  > ### Title: Distribution of LR Test Statistic for the Rank Test
  > ### Aliases: FCVARbootRank
  >
  > ### ** Examples
  >
  > opt <- FCVARoptions()
> opt$gridSearch   <- 0 # Disable grid search in optimization.
> opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
> opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
> opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
> DefaultOpt$plotRoots <- 0
Error in DefaultOpt$plotRoots <- 0 : object 'DefaultOpt' not found
Execution halted
√  checking for unstated dependencies in 'tests' ...
-  checking tests ...
√  Running 'testthat.R' [45s] (45.5s)

See
'C:/Users/le279259/AppData/Local/Temp/1/RtmpeWh3xu/FCVAR.Rcheck/00check.log'
for details.

-- R CMD check results ------------------------------------- FCVAR 0.0.0.9000 ----
  Duration: 1m 33.4s

> checking examples ... ERROR
Running examples in 'FCVAR-Ex.R' failed
The error most likely occurred in:

  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: FCVARbootRank
  > ### Title: Distribution of LR Test Statistic for the Rank Test
  > ### Aliases: FCVARbootRank
  >
  > ### ** Examples
  >
  > opt <- FCVARoptions()
> opt$gridSearch   <- 0 # Disable grid search in optimization.
> opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
> opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
> opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
> DefaultOpt$plotRoots <- 0
Error in DefaultOpt$plotRoots <- 0 : object 'DefaultOpt' not found
Execution halted

> checking whether package 'FCVAR' can be installed ... WARNING
See below...

> checking index information ... WARNING
Demos with missing or empty index information:
  FCVAR_extension_JNP2014
FCVAR_replication_JNP2014
See sections 'The INDEX file' and 'Package subdirectories' in the
'Writing R Extensions' manual.

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

    See section 'Generic functions and methods' in the 'Writing R
  Extensions' manual.

> checking Rd \usage sections ... WARNING
Undocumented arguments in documentation object 'FCVARlikeFull'
'beta' 'rho'
Documented arguments not in \usage in documentation object 'FCVARlikeFull':
  'betaHat' 'rhoHat'

Undocumented arguments in documentation object 'GetFreeParams'
'p'
Documented arguments not in \usage in documentation object 'GetFreeParams':
  'x'

Functions with \usage entries need to have the appropriate \alias
entries, and all their arguments documented.
The \usage entries must correspond to syntactically valid R code.
See chapter 'Writing R documentation files' in the 'Writing R
  Extensions' manual.

> checking data for ASCII and uncompressed saves ... WARNING
Warning: package needs dependence on R (>= 2.10)

> checking top-level files ... NOTE
Non-standard file/directory found at top level:
  'article'

> checking R code for possible problems ... NOTE
FCVARestn: no visible global function definition for 'optim'
FCVARestn: no visible global function definition for 'constrOptim'
FCVARhypoTest: no visible global function definition for 'pchisq'
FCVARlagSelect: no visible global function definition for 'pchisq'
FCVARlikeGrid: no visible global function definition for 'optim'
FCVARoptionUpdates: no visible binding for global variable 'end'
FCVARsim: no visible global function definition for 'rnorm'
FCVARsimBS: no visible global function definition for 'rnorm'
FCVARsimBS: possible error in matrix(1, nrow = T, nol = 1): unused
argument (nol = 1)
FracDiff: no visible global function definition for 'nextn'
FracDiff: no visible global function definition for 'fft'
LMtest: no visible global function definition for 'pchisq'
Qtest: no visible global function definition for 'pchisq'
plot.FCVARlikeGrid: no visible global function definition for 'pdf'
plot.FCVARlikeGrid: no visible global function definition for 'png'
plot.FCVARlikeGrid: no visible global function definition for 'rainbow'
plot.FCVARlikeGrid: no visible global function definition for 'par'
plot.FCVARlikeGrid: no visible global function definition for 'persp'
plot.FCVARlikeGrid: no visible global function definition for 'plot'
plot.FCVARlikeGrid: no visible global function definition for 'dev.off'
plot.GetCharPolyRoots: no visible global function definition for 'pdf'
plot.GetCharPolyRoots: no visible global function definition for 'png'
plot.GetCharPolyRoots: no visible global function definition for 'plot'
plot.GetCharPolyRoots: no visible global function definition for
'lines'
plot.GetCharPolyRoots: no visible global function definition for
'points'
plot.GetCharPolyRoots: no visible global function definition for
'dev.off'
print.MVWNtest: no visible binding for global variable 'p'
Undefined global functions or variables:
  constrOptim dev.off end fft lines nextn optim p par pchisq pdf persp
plot png points rainbow rnorm
Consider adding
importFrom("grDevices", "dev.off", "pdf", "png", "rainbow")
importFrom("graphics", "lines", "par", "persp", "plot", "points")
importFrom("stats", "constrOptim", "end", "fft", "nextn", "optim",
           "pchisq", "rnorm")
to your NAMESPACE file.

1 error x | 5 warnings x | 2 notes x
>
