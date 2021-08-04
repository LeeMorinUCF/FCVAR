# FCVAR
An R Package for the Fractionally Cointegrated VAR Model

## Description

Estimation and inference using the Fractionally Cointegrated 
Vector Autoregressive (VAR) model. It includes functions for model specification, 
including lag selection and cointegration rank selection, as well as a comprehensive
set of options for hypothesis testing, including tests of hypotheses on the 
cointegrating relations, the adjustment coefficients and the fractional 
differencing parameters. 
See the file ```FCVAR_README.pdf``` for examples
and the Webpage https://sites.google.com/view/mortennielsen/software
for more information about the FCVAR model.


## How to Install FCVAR

Install the development version on GitHub using the ```devtools``` package:

```
library(devtools)
devtools::install_github("LeeMorinUCF/FCVAR")
```



## Release Notes

# August 4, 2021: Fourth submission.

* Added link to reference paper in the DESCRIPTION file. 
* Changed examples with ```\dontrun``` to ```\donttest``` for examples
with run time than took longer than 5s.
* Removed example from demo that changed ```par()``` settings.
* For function ```plot.FCVAR_grid()``` that changes ```par()``` settings, 
because it creates a figure with thinner margins, 
inserted command ```on.exit(par(oldpar))``` to restore user's settings, 
immediately after the change to ```par()```. 


# July 30, 2021: Third submission.

* Excluded examples with ```\dontrun``` if run time took longer than 5s.


# July 30, 2021: Second submission.

* Added cran-comments.md to .Rbuildignore.



### July 30, 2021,  ```FCVAR_0.1.0```: 

* First submission of the FCVAR package. 
* Added a `NEWS.md` file to track changes to the package.


#### Checks on Rhub:
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit

#### Checks on local machine:
* Windows 10 Enterprise, version 20H2, OS build 19042.985, R 4.0.5

#### Checks on win-builder:
* devel and release


#### R CMD check results

There were no ERRORs or WARNINGs.

There were three NOTEs:

* New submission

* Authors' names, conjugations of the verb "cointegrate," the verb "difference," 
and the acronym "FCVAR" are not mis-spelled. 

* Elapsed time > 10s for some examples.


#### Downstream dependencies

There are currently no downstream dependencies for this package. 

