# FCVAR 0.1.3


### May 4, 2022: ```FCVAR v0.1.3```

* Adjustment to test case to prevent false alarm 
from tests run on M1mac platform. 


## Test environments


### Checks on Rhub:
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Windows Server 2022, R-devel, 64 bit

### Checks on local machine:
* Windows 10 Enterprise, version 20H2, OS build 19042.985, R 4.0.5

### Checks on win-builder:
* devel and release

### Checks on M1mac
* using R version 4.2.0 RC (2022-04-15 r82193)
* using platform: aarch64-apple-darwin20 (64-bit)



## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE from the check on Windows Server 2022, R-devel, 64 bit on Rhub:

* checking for detritus in the temp directory ... NOTE
  'lastMiKTeXException'
* Apparently, this is a known issue with Rhub and does not suggest a problem with the package. 


## Downstream dependencies

There are currently no downstream dependencies for this package. 




## Previous Submissions


### March 10, 2022: ```FCVAR v0.1.2```

* Adjustment for issue with counting the number of free parameters in constrained model. 
* Added warning message for rank test with rank above 12. In this case, 
p-values are only available by simulation with the ```FCVARbootRank()``` function. 


### August 6, 2021: ```FCVAR v0.1.1```

* For models with restricted estimation, skipped tests on CRAN 
that compare printed output and require exact equality (i.e. without tolerance), 
but retained tests on same with default tolerance for numerical equality. 


### August 5, 2021: ```FCVAR v0.1.1```

* Modified tolerance of bootstrap test results for Solaris platforms.
Test results are within acceptable tolerance. 

### August 4, 2021: Fourth submission (```FCVAR v0.1.0```).

* Added link to reference paper in the DESCRIPTION file. 
* Changed examples with ```\dontrun``` to ```\donttest``` for examples
with run time than took longer than 5s.
* Removed example from demo that changed ```par()``` settings.
* For function ```plot.FCVAR_grid()``` that changes ```par()``` settings, 
because it creates a figure with thinner margins, 
inserted command ```on.exit(par(oldpar))``` to restore user's settings, 
immediately after the change to ```par()```. 

### July 30, 2021: Third submission (```FCVAR v0.1.0```).

* Excluded examples with ```\dontrun``` if run time took longer than 5s.


### July 30, 2021: Second submission (```FCVAR v0.1.0```).

* Added cran-comments.md to .Rbuildignore.


### July 30, 2021: First submission (```FCVAR v0.1.0```).

First submission of the FCVAR package.  

