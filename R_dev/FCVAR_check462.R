
# Results from check number 462:

-- R CMD check results ------------------------------------------------------------ FCVAR 0.0.0.9000 ----
  Duration: 4m 36.3s

0 errors √ | 0 warnings √ | 0 notes √


# Output from devtools function calls:

> devtools::load_all()
i Loading FCVAR

# Trimmed examples to avoid skipping tests.

> devtools::test()
i Loading FCVAR
i Testing FCVAR
√ |  OK F W S | Context
√ |  12       | Estimation [33.3 s]
√ |   9       | Postestimation [13.3 s]
√ |   6       | Specification [15.6 s]

== Results ==============================================================================================
  Duration: 62.2 s

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 27 ]

# All pass and complete within about a minute.

> devtools::check()
i Updating FCVAR documentation
i Loading FCVAR
Writing NAMESPACE
Writing NAMESPACE
Writing plot.FCVAR_grid.Rd
-- Building ------------------------------------------------------------------------------------ FCVAR --
  Setting env vars:
  * CFLAGS    : -Wall -pedantic
* CXXFLAGS  : -Wall -pedantic
* CXX11FLAGS: -Wall -pedantic
---------------------------------------------------------------------------------------------------------
  √  checking for file 'C:\Users\le279259\OneDrive - University of Central Florida\Documents\Research\FCVAR\GitRepo\FCVAR/DESCRIPTION'
-  preparing 'FCVAR': (4.4s)
√  checking DESCRIPTION meta-information ...
-  checking for LF line-endings in source and make files and shell scripts (2s)
-  checking for empty or unneeded directories
Removed empty directory 'FCVAR/vignettes'
-  building 'FCVAR_0.0.0.9000.tar.gz'

-- Checking ------------------------------------------------------------------------------------ FCVAR --
  Setting env vars:
  * _R_CHECK_CRAN_INCOMING_REMOTE_: FALSE
* _R_CHECK_CRAN_INCOMING_       : FALSE
* _R_CHECK_FORCE_SUGGESTS_      : FALSE
* NOT_CRAN                      : true
-- R CMD check ------------------------------------------------------------------------------------------
  -  using log directory 'C:/Users/le279259/AppData/Local/Temp/1/RtmpYb9ALL/FCVAR.Rcheck'
-  using R version 4.0.5 (2021-03-31)
-  using platform: x86_64-w64-mingw32 (64-bit)
-  using session charset: ISO8859-1
-  using options '--no-manual --as-cran'
√  checking for file 'FCVAR/DESCRIPTION'
-  this is package 'FCVAR' version '0.0.0.9000'
-  package encoding: UTF-8
√  checking package namespace information ...
√  checking package dependencies (1s)
√  checking if this is a source package ...
√  checking if there is a namespace
√  checking for .dll and .exe files
√  checking for hidden files and directories ...
√  checking for portable file names ...
√  checking whether package 'FCVAR' can be installed (4.9s)
√  checking package directory
√  checking for future file timestamps (377ms)
√  checking DESCRIPTION meta-information (347ms)
√  checking top-level files
√  checking for left-over files
√  checking index information ...
√  checking package subdirectories ...
√  checking R files for non-ASCII characters ...
√  checking R files for syntax errors ...
√  checking whether the package can be loaded ...
√  checking whether the package can be loaded with stated dependencies ...
√  checking whether the package can be unloaded cleanly ...
√  checking whether the namespace can be loaded with stated dependencies ...
√  checking whether the namespace can be unloaded cleanly ...
√  checking loading without being on the library search path (404ms)
√  checking dependencies in R code ...
√  checking S3 generic/method consistency (552ms)
√  checking replacement functions ...
√  checking foreign function calls ...
√  checking R code for possible problems (6.5s)
√  checking Rd files (437ms)
√  checking Rd metadata ...
√  checking Rd line widths ...
√  checking Rd cross-references ...
√  checking for missing documentation entries ...
√  checking for code/documentation mismatches (659ms)
√  checking Rd \usage sections (774ms)
√  checking Rd contents ...
√  checking for unstated dependencies in examples (443ms)
√  checking contents of 'data' directory ...
√  checking data for non-ASCII characters ...
√  checking data for ASCII and uncompressed saves ...
√  checking examples (3m 15.9s)
Examples with CPU (user + system) or elapsed time > 5s
user system elapsed
FCVARhypoTest       30.53   0.06   30.60
FCVARestn           30.13   0.00   30.14
plot.FCVAR_grid     13.00   0.00   13.02
FCVARlikeGrid       12.75   0.02   12.80
FCVARboot           12.66   0.01   12.68
FCVARbootRank        9.73   0.00    9.73
FCVARforecast        8.13   0.00    8.13
FCVARsim             7.68   0.00    7.69
FCVARsimBS           7.63   0.02    7.64
MVWNtest             7.39   0.00    7.39
plot.FCVAR_roots     7.39   0.00    7.39
summary.MVWN_stats   7.38   0.00    7.37
GetCharPolyRoots     7.15   0.00    7.19
summary.FCVAR_model  7.14   0.01    7.16
summary.FCVAR_roots  7.06   0.00    7.07
FCVARrankTests       5.59   0.01    5.62
summary.FCVAR_ranks  5.19   0.00    5.18
√  checking for unstated dependencies in 'tests' ...
-  checking tests ...
√  Running 'testthat.R' [58s] (58.1s)
√  checking for non-standard things in the check directory (58.1s)
√  checking for detritus in the temp directory


-- R CMD check results ------------------------------------------------------------ FCVAR 0.0.0.9000 ----
  Duration: 4m 36.3s

0 errors √ | 0 warnings √ | 0 notes √


# All pass.


# Pre-submission procedures:

# Sequence of checks for submission to CRAN

# Rhub checks on Windows and Two Linux platforms:
# rhub_results <- rhub::check_for_cran()
# rhub_results$cran_summary()
# devtools::check_win_release()
# devtools::check_win_devel()

# When all checks pass, update cran-comments.md,
# NEWS.md, and README.md.

# Then, submit to CRAN:

# This way, you get a series of questions to
# run through a checklist:
# devtools::release()

# This way is a short cut,
# which skips all the questions:
# devtools::submit_cran()


# Then wait anxiously for a response.





# Now make more rigorous tests for submission.


# > library(rhub)
#
# > rhub_results <- rhub::check_for_cran()
#
# rhub_results$cran_summary()


# Found problems related to differences in numerical precision across platforms.
# expect_equal(10, 10)
# expect_equal(10, 10 + 1e-7)
# expect_equal(10, 10 + 1e-6)
# expect_equal(10, 11)



# Common sequence of builds, tests and checks.
devtools::load_all()
devtools::document()
devtools::build_manual()
# Run these locally first to avoid embarrassment.
devtools::test()
devtools::check()
# Run these first, since they run quickly but return results with a delay.
devtools::check_win_release()
devtools::check_win_devel()
# Run this next, since it takes longer but return results quickly.
rhub_results <- rhub::check_for_cran()
devtools::release()
# devtools::submit_cran() # Shortcut without all the questions.


# After the release, an automated file CRAN-RELEASES will appear in the main directory.
# Create a tag in git as follows:
# git tag name_of_tag commit_hash
# Then push the tag, like any other commit, with this command:
# git push origin name_of_tag
# name_of_tag is usually of the form package_name_v9.9.9

# Then go to github.com and create a release linked to this tag.
# Copy the description from your NEWS.md file in the release description.

# Then delete the CRAN-RELEASES file, which will be replaced on the next release.



