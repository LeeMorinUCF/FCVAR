

# Results from check 3:

nickname       Feather Spray
> devtools::check()
Updating FCVAR documentation
Writing NAMESPACE
Loading FCVAR
Writing NAMESPACE
Writing FCVARestn.Rd
Writing print.MVWNtest.Rd
-- Building ----------------------------------------------------------- FCVAR --
  Setting env vars:
  * CFLAGS    : -Wall -pedantic
* CXXFLAGS  : -Wall -pedantic
* CXX11FLAGS: -Wall -pedantic
--------------------------------------------------------------------------------
  √  checking for file 'C:\Users\le279259\Documents\Research\FCVAR\GitRepo\FCVAR/DESCRIPTION' ...
-  preparing 'FCVAR': (12.8s)
√  checking DESCRIPTION meta-information ...
-  checking whether 'INDEX' is up-to-date ... NO (591ms)
-  use '--force' to remove the existing 'INDEX'
-  checking for LF line-endings in source and make files and shell scripts (1.4s)
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
-  using options '--no-manual --as-cran'
√  checking for file 'FCVAR/DESCRIPTION' ...
-  this is package 'FCVAR' version '0.0.0.9000'
-  package encoding: UTF-8
√  checking package namespace information ...
√  checking package dependencies (989ms)
√  checking if this is a source package ...
√  checking if there is a namespace
√  checking for executable files (416ms)
√  checking for hidden files and directories ...
√  checking for portable file names ...
√  checking whether package 'FCVAR' can be installed (4.3s)
√  checking installed package size ...
√  checking package directory
W  checking DESCRIPTION meta-information (340ms)
Dependence on R version '3.5.1' not with patchlevel 0
√  checking top-level files
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
W  checking S3 generic/method consistency (512ms)
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
√  checking R code for possible problems (5.6s)
√  checking Rd files (513ms)
√  checking Rd metadata ...
√  checking Rd line widths ...
√  checking Rd cross-references ...
√  checking for missing documentation entries ...
√  checking for code/documentation mismatches (730ms)
√  checking Rd \usage sections (824ms)
√  checking Rd contents ...
√  checking for unstated dependencies in examples (614ms)
√  checking contents of 'data' directory
√  checking data for non-ASCII characters ...
√  checking data for ASCII and uncompressed saves ...
E  checking examples (52.6s)
Running examples in 'FCVAR-Ex.R' failed
The error most likely occurred in:

  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: FCVARestn
  > ### Title: Estimate FCVAR model
  > ### Aliases: FCVARestn
  >
  > ### ** Examples
  >
  > opt <- FCVARoptions()
> opt$gridSearch   <- 0 # Disable grid search in optimization.
> opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
> opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
> opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
> x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
> m1 <- FCVARestn(x, k = 2, r = 1, opt)

--------------------------------------------------------------------------------
  Fractionally Cointegrated VAR: Estimation Results
--------------------------------------------------------------------------------
  Dimension of system:       3      Number of observations in sample:          316
Number of lags:            2      Number of observations for estimation:     316
Restricted constant:      No      Initial values:                              0
Unrestricted constant:    No      Level parameter:                           Yes
Starting value for d:    1.000    Parameter space for d: (0.010 , 2.000)
Starting value for b:    1.000    Parameter space for b: (0.010 , 2.000)
--------------------------------------------------------------------------------
  Cointegrating rank:            1  AIC:              -848.348
Log-likelihood:          451.174  BIC:              -746.943
log(det(Omega_hat)):     -11.369  Free parameters:        27
--------------------------------------------------------------------------------
  Fractional parameters:
  --------------------------------------------------------------------------------
  Coefficient               Estimate                Standard error
--------------------------------------------------------------------------------
  d                       0.569                      0.049
--------------------------------------------------------------------------------
  --------------------------------------------------------------------------------
  Cointegrating equations (beta):
  --------------------------------------------------------------------------------
  Variable        CI equation 1
--------------------------------------------------------------------------------
  Var1              1.000
Var2              0.111
Var3             -0.240
--------------------------------------------------------------------------------
  Note: Identifying restriction imposed.
--------------------------------------------------------------------------------
  Adjustment matrix (alpha):
  --------------------------------------------------------------------------------
  Variable        CI equation 1
--------------------------------------------------------------------------------
  Var 1            -0.180
SE 1         (   0.064  )
Var 2             0.167
SE 2         (   0.194  )
Var 3             0.037
SE 3         (   0.014  )
--------------------------------------------------------------------------------
  Note: Standard errors in parenthesis.
--------------------------------------------------------------------------------
  Long-run matrix (Pi):
  --------------------------------------------------------------------------------
  Variable         Var 1          Var 2          Var 3
--------------------------------------------------------------------------------
  Var 1           -0.180         -0.020          0.043
Var 2            0.167          0.019         -0.040
Var 3            0.037          0.004         -0.009
--------------------------------------------------------------------------------

  --------------------------------------------------------------------------------
  Level parameter (mu):
  --------------------------------------------------------------------------------
  Var 1            -0.345
SE 1         (   0.069  )
Var 2            11.481
SE 2         (   0.548  )
Var 3            -2.872
SE 3         (   0.033  )
--------------------------------------------------------------------------------
  Note: Standard errors in parenthesis (from numerical Hessian)
but asymptotic distribution is unknown.
--------------------------------------------------------------------------------
  Lag matrix 1 (Gamma_1):
  --------------------------------------------------------------------------------
  Variable         Var 1          Var 2          Var 3
--------------------------------------------------------------------------------
  Var 1            0.276         -0.032         -0.510
SE 1        (   0.160  )   (   0.026  )   (   0.513  )
Var 2           -0.148          1.126         -3.285
SE 2        (   0.378  )   (   0.196  )   (   1.975  )
Var 3           -0.052          0.008          0.711
SE 3        (   0.022  )   (   0.005  )   (   0.170  )
--------------------------------------------------------------------------------
  Note: Standard errors in parentheses.
--------------------------------------------------------------------------------
  Lag matrix 2 (Gamma_2):
  --------------------------------------------------------------------------------
  Variable         Var 1          Var 2          Var 3
--------------------------------------------------------------------------------
  Var 1            0.566          0.106          0.609
SE 1        (   0.182  )   (   0.045  )   (   0.612  )
Var 2            0.493         -0.462          0.450
SE 2        (   0.562  )   (   0.198  )   (   2.627  )
Var 3           -0.039         -0.020          0.318
SE 3        (   0.032  )   (   0.008  )   (   0.143  )
--------------------------------------------------------------------------------
  Note: Standard errors in parentheses.
--------------------------------------------------------------------------------
  --------------------------------------------------------------------------------
  Roots of the characteristic polynomial
--------------------------------------------------------------------------------
  Number     Real part    Imaginary part       Modulus
--------------------------------------------------------------------------------
  1         -2.893         -0.000            2.893
2         -1.522         -0.000            1.522
3          1.010         -0.927            1.371
4          1.010          0.927            1.371
5          1.108          0.000            1.108
6          1.000          0.000            1.000
7          1.000          0.000            1.000
8          0.944         -0.261            0.980
9          0.944          0.261            0.980
--------------------------------------------------------------------------------

  --------------------------------------------------------------------------------
  Restrictions imposed on the following parameters:
  - Psi. For details see "options$R_psi"
--------------------------------------------------------------------------------

  >
  > opt1 <- opt
> opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
> opt1$r_psi <- 1
> m1r1 <- FCVARestn(x, k, r, opt1)
Error in Lbk(x, b, k) : object 'k' not found
Calls: FCVARestn ... FCVARlikeMu -> GetParams -> TransformData -> FracDiff -> Lbk
Execution halted
√  checking for unstated dependencies in 'tests' ...
-  checking tests ...
√  Running 'testthat.R' [49s] (48.7s)

See
'C:/Users/le279259/AppData/Local/Temp/1/RtmpeWh3xu/FCVAR.Rcheck/00check.log'
for details.

-- R CMD check results ------------------------------------- FCVAR 0.0.0.9000 ----
  Duration: 2m 2s

> checking examples ... ERROR
Running examples in 'FCVAR-Ex.R' failed
The error most likely occurred in:

  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: FCVARestn
  > ### Title: Estimate FCVAR model
  > ### Aliases: FCVARestn
  >
  > ### ** Examples
  >
  > opt <- FCVARoptions()
> opt$gridSearch   <- 0 # Disable grid search in optimization.
> opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
> opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
> opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
> x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
> m1 <- FCVARestn(x, k = 2, r = 1, opt)

--------------------------------------------------------------------------------
  Fractionally Cointegrated VAR: Estimation Results
--------------------------------------------------------------------------------
  Dimension of system:       3      Number of observations in sample:          316
Number of lags:            2      Number of observations for estimation:     316
Restricted constant:      No      Initial values:                              0
Unrestricted constant:    No      Level parameter:                           Yes
Starting value for d:    1.000    Parameter space for d: (0.010 , 2.000)
Starting value for b:    1.000    Parameter space for b: (0.010 , 2.000)
--------------------------------------------------------------------------------
  Cointegrating rank:            1  AIC:              -848.348
Log-likelihood:          451.174  BIC:              -746.943
log(det(Omega_hat)):     -11.369  Free parameters:        27
--------------------------------------------------------------------------------
  Fractional parameters:
  --------------------------------------------------------------------------------
  Coefficient               Estimate                Standard error
--------------------------------------------------------------------------------
  d                       0.569                      0.049
--------------------------------------------------------------------------------
  --------------------------------------------------------------------------------
  Cointegrating equations (beta):
  --------------------------------------------------------------------------------
  Variable        CI equation 1
--------------------------------------------------------------------------------
  Var1              1.000
Var2              0.111
Var3             -0.240
--------------------------------------------------------------------------------
  Note: Identifying restriction imposed.
--------------------------------------------------------------------------------
  Adjustment matrix (alpha):
  --------------------------------------------------------------------------------
  Variable        CI equation 1
--------------------------------------------------------------------------------
  Var 1            -0.180
SE 1         (   0.064  )
Var 2             0.167
SE 2         (   0.194  )
Var 3             0.037
SE 3         (   0.014  )
--------------------------------------------------------------------------------
  Note: Standard errors in parenthesis.
--------------------------------------------------------------------------------
  Long-run matrix (Pi):
  --------------------------------------------------------------------------------
  Variable         Var 1          Var 2          Var 3
--------------------------------------------------------------------------------
  Var 1           -0.180         -0.020          0.043
Var 2            0.167          0.019         -0.040
Var 3            0.037          0.004         -0.009
--------------------------------------------------------------------------------

  --------------------------------------------------------------------------------
  Level parameter (mu):
  --------------------------------------------------------------------------------
  Var 1            -0.345
SE 1         (   0.069  )
Var 2            11.481
SE 2         (   0.548  )
Var 3            -2.872
SE 3         (   0.033  )
--------------------------------------------------------------------------------
  Note: Standard errors in parenthesis (from numerical Hessian)
but asymptotic distribution is unknown.
--------------------------------------------------------------------------------
  Lag matrix 1 (Gamma_1):
  --------------------------------------------------------------------------------
  Variable         Var 1          Var 2          Var 3
--------------------------------------------------------------------------------
  Var 1            0.276         -0.032         -0.510
SE 1        (   0.160  )   (   0.026  )   (   0.513  )
Var 2           -0.148          1.126         -3.285
SE 2        (   0.378  )   (   0.196  )   (   1.975  )
Var 3           -0.052          0.008          0.711
SE 3        (   0.022  )   (   0.005  )   (   0.170  )
--------------------------------------------------------------------------------
  Note: Standard errors in parentheses.
--------------------------------------------------------------------------------
  Lag matrix 2 (Gamma_2):
  --------------------------------------------------------------------------------
  Variable         Var 1          Var 2          Var 3
--------------------------------------------------------------------------------
  Var 1            0.566          0.106          0.609
SE 1        (   0.182  )   (   0.045  )   (   0.612  )
Var 2            0.493         -0.462          0.450
SE 2        (   0.562  )   (   0.198  )   (   2.627  )
Var 3           -0.039         -0.020          0.318
SE 3        (   0.032  )   (   0.008  )   (   0.143  )
--------------------------------------------------------------------------------
  Note: Standard errors in parentheses.
--------------------------------------------------------------------------------
  --------------------------------------------------------------------------------
  Roots of the characteristic polynomial
--------------------------------------------------------------------------------
  Number     Real part    Imaginary part       Modulus
--------------------------------------------------------------------------------
  1         -2.893         -0.000            2.893
2         -1.522         -0.000            1.522
3          1.010         -0.927            1.371
4          1.010          0.927            1.371
5          1.108          0.000            1.108
6          1.000          0.000            1.000
7          1.000          0.000            1.000
8          0.944         -0.261            0.980
9          0.944          0.261            0.980
--------------------------------------------------------------------------------

  --------------------------------------------------------------------------------
  Restrictions imposed on the following parameters:
  - Psi. For details see "options$R_psi"
--------------------------------------------------------------------------------

  >
  > opt1 <- opt
> opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
> opt1$r_psi <- 1
> m1r1 <- FCVARestn(x, k, r, opt1)
Error in Lbk(x, b, k) : object 'k' not found
Calls: FCVARestn ... FCVARlikeMu -> GetParams -> TransformData -> FracDiff -> Lbk
Execution halted

> checking DESCRIPTION meta-information ... WARNING
Dependence on R version '3.5.1' not with patchlevel 0

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

1 error x | 3 warnings x | 0 notes √
>
