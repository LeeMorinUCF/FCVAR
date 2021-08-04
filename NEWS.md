# FCVAR 0.1.0

# August 4, 2021: Fourth submission.

* Added link to reference paper in the DESCRIPTION file. 
* Changed examples with ```\dontrun``` to ```\donttest``` for examples
with run time than took longer than 5s.
* Removed example from demo that changed ```par()```.
* For function ```plot.FCVAR_grid()``` that changes ```par()``` settings, 
because it creates a figure with thinner margins, 
inserted command ```on.exit(par(oldpar))``` to restore user's settings. 


# July 30, 2021: Third submission.

* Excluded examples with ```\dontrun``` if run time took longer than 5s.


# July 30, 2021: Second submission.

* Added cran-comments.md to .Rbuildignore.


## July 30, 2021: 

First submission of the FCVAR package. 
