/* 
Program: Fast Fractional Differencing
Author:  Michal K. Popiel (www.michalpopiel.com)
Version: 1.0 
Date: 	 April 5, 2014 
Description: uses a fast fractional differencing algorithm developed by Jensen and Nielsen 
	     (2013) to fractionally difference a series
*/

program fracdiff
	version 12
	syntax varname(ts max =1) [if] [in], D(real) Name(string)
	
	quietly gen `name' = .
	marksample touse
	mata: fd("`varlist'", "`touse'", `d', "`name'")

end
version 12
mata:

/* This function computes the cumulative product of the columns in a matrix */
real matrix cumprod(real matrix x){
	nrow = rows(x)
	ncol = cols(x)
	A = x
	for (i=1; i<=ncol; i++){
		for(j=2; j<=nrow; j++){
			A[j,i] = A[j-1,i]*x[j,i]
			}
		}
	return (A)
}

/* Thus function tiles a matrix "x" p times horizontally and n times vertically */
real matrix repmat(real matrix A, p, n){
	x=A
	for (i=1; i<=n; i++){
		x = x,x
		}
	for(j=1; j<=p; j++){
		x = x\x
		}
	return (x)
}



// this function implements a fast fractional algorithm to fractionally difference a series stored
// in matrix x, the series are the rows, the d is the order of integration.
real matrix fastfd(real matrix x, real d){
	
	nrow = rows(x)
	ncol = cols(x)
	
	k = (1..nrow-1)'
	
	b = (k:-d:-1):/k
	b = (1)\cumprod(b)\J(nrow,1,0)
	ftb = fft( ftpad(b) )
	fdx = x
	
	for(i=1; i<=ncol; i++){
		dx = invfft( ftb:*fft( ftpad(x[.,i]\J(nrow,1,0)) ) )
		fdx[.,i] = Re(dx[1..nrow,.])
	}
	
	return(fdx)
}

void fd(string scalar varname, string scalar touse, real d, string scalar fdvarname)
{
	x = .;
	// call variable to fractionally difference
	st_view(x, ., varname, touse)
	// perform the calculation
	y = fastfd(x,d);
	// output the result to a new variable specified by the user
	st_store(.,fdvarname, y)
}

end


