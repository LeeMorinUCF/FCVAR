

/*
  This is an example do-file for using "FCVAR.ado" and "fracdiff.ado"
  Author:  Michal K. Popiel (www.michalpopiel.com) 
  Date:    April 7, 2014 
*/

clear


/*
   Import data from http://people.exeter.ac.uk/jehd201/bdpdata.txt
   This data is on popularity of the two major political parties in England
   Reference:
   Byers, D., Davidson, J., and Peel, D. (1997). Modelling political popularity: an analysis
	of long-range dependence in opinion poll series. Journal of the Royal Statistical Society:
	Series A (Statistics in Society), 160(3):471-490.
*/
use "bdp97_data.dta"

// look at individual series to see if they are of close fractional order
FCVAR logl, l(0) r(0) c(0)
return list

// generate the fractionally differenced series and plot it
fracdiff logl, d(`r(d)') n(fdlogl)

tsline fdlogl

FCVAR logc, l(0) r(0) c(0)

// generate the fractionally differenced series and plot it
fracdiff logc, d(`r(d)') n(fdlogc)

tsline fdlogc


// look at fractionally cointegrated system

// likelihood ratio test for lag selection (example)

// Unrestricted Model:
FCVAR logc logl, l(2) r(1) c(1)

scalar logUNR = r(loglike)
// Restricted Model:
FCVAR logc logl, l(1) r(1) c(1)


scalar logR = r(loglike)

// Test:
scalar LRstat = 2*(logUNR - logR)
di "chi2(4) = " LRstat
di "Prob > chi2 = "chi2tail(4, LRstat)

// look at residuals
* matrix list  r(resid)

// store the residuals 
matrix epsilon = r(resid)

svmat epsilon

// plot the residuals
tsline epsilon1 epsilon2 

