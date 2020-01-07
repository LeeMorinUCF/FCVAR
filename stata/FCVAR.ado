/* 
Program: Fractionally Cointegrated Vector Autoregressive Model Estimation
Author:  Michal K. Popiel (www.michalpopiel.com)
Version: 1.0 
Date: 	 April 5, 2014 
* please cite when using results produced by this program *
Description: This program estimates the FCVAR model developed by Nielsen and Johansen (Econometrica 2012)
- Currently the model does not support initial values
- Only the case d=b can be estimated
- Restriction testing has not yet been implemented
- Maximization of the likelihood occurs in two-steps to increase accuracy (see comments below)
- The program closely follows the Matlab package developed by Nielsen and Lee (2013)
*/


program FCVAR, rclass
	version 12
	syntax varlist(ts) [if] [in], Lags(real) Rank(real) Constant(real) ///
			[dmin(real 0.1) dmax(real 1.5) dstep(real 0.01)]
	
	* create some variables to be returned after estimation
	global matrix estimates
	global matrix residuals
	global matrix se
	global matrix logl
	global scalar dstar
	
	
	
	marksample touse
	
	* run the main estimation program
	mata: fcvarRun("`varlist'", "`touse'", `lags', `rank', `constant' ///
			, `dmin', `dmax', `dstep')


	return matrix est = estimates
	return matrix resid = residuals
	return matrix stderror = se
	return scalar d = dstar
	return scalar loglike = logl[1,1]
end


// The following is the MATA code that runs the estimation
version 12
mata:

/* this function serves as an interface with the ado file and mata functions
   it takes the inputs specified by the user, estimates and returns all of the
   coefficient estimates
*/ 
void fcvarRun(varname, touse, real k, real r, real m, real dmin, real dmax, real dstep){

	x = .;

	st_view(x, ., varname, touse)
	
	// ESTIMATE THE MODEL
	d = SolveModelGrid(x,k,r,m,dmin,dmax,dstep);
	
	st_numscalar("dstar", d);

	
	// RECOVER THE PARAMETERS	
	est = GetParams(x,k,r,d,m);
	est = Re (est);
	st_matrix("estimates", est);
	
	p = cols(x);
	nCols = cols(est)
	
	i = r
	if(r==0){
		i = 1
	}

	alphaHat = est[1..p, (p+1)..(p+1)+i-1];
	betaHat  = est[1..p, (p+1+i)..(p+1+i)+i-1];
	rhoHat   = est[1, (p+1+2*i)..(p+1+2*i)+i-1];
	gammaHat = est[1..p, (p+1+3*i)..nCols];
	
	// OBTAIN THE RESIDUALS
	res = GetResiduals(x, k, r, d, m, alphaHat, betaHat, rhoHat, gammaHat);
	
	st_matrix("residuals", res);

	// CALCULATE STANDARD ERRORS	
	H = FCVARhess(x, k, r, d, m, est);
	SE = sqrt(-diagonal(luinv(H)));
	
	st_matrix("se", SE);

}

/* This function computes the cumulative product of the columns in a matrix 
   It is used in by the fractional differencing algorithm
*/
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

/* Thus function tiles a matrix "x" p times horizontally and n times vertically 
   It is meant to be similar to the repmat matrix in matlab
*/
real matrix repmat(real matrix A, p, n){
	x=A
	for (i=1; i<n; i++){
		x = x,A
		}
	for(j=1; j<p; j++){
		x = x\A
		}
	return (x)
}



/* this function implements a fast fractional algorithm to fractionally difference a series stored
   in matrix x, the series are the rows, the d is the order of integration.
   the algorithm is based on Nielsen and Jensen (2013)
   see the paper for details
*/ 

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

/*
  This function can be used stand alone to interface with an ado file and return
  a fractionally differenced series to a new variable specified by the user
*/
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

/*
  This function Lbk returns the lagged polynomial of a matrix
  It stackes the lagged matrices next to each other
*/
real matrix lbk(real matrix x, real d, k){
	// number of variables
	p = cols(x) 
	// initialize the return matrix
	lbkx = repmat(x,0,k) 
	if(k>0){
		lbkx[.,1..p] = lbkx[.,1..p] - fastfd(lbkx[.,1..p],d)
		for(i=2; i<=k; i++){
			   lbkx[.,(p*(i-1)+1)..(p*i)] = lbkx[.,(p*(i-2)+1)..(p*(i-1))] - fastfd(lbkx[.,(p*(i-2)+1)..(p*(i-1))],d)
			}
		}
	return(lbkx)
}


/* 
   This function transforms the data so that the reduced rank regression can
   be run more efficiently. Z0 is the fractionally differenced series. Z1 is the 
   fractionally lagged series with a column of ones if a constant is included in the 
   user-specified entry. Z2 is the fractionally differenced matrix of fractional lags.
   To improve speed and memory usage, pointers are used.
*/
void TransformData(x, k, d, m, pointer Z0, Z1, Z2){
	nrow = rows(x);
	*Z0 = fastfd(x,d);
	*Z2 = fastfd(lbk(x,d,k),d);
	if(m==1){
	  A = x,J(nrow,1,1);
	}
	else{
	  A = x;
	}
	*Z1 = lbk(A,d,1);
}


/*
   This function concentrates out the parameters and performs the reduced rank regression to 
   solve for the OmegaHat matrix that enters the likelihood estimation. It returns
   the estimated parameters for a given x,k,r,d,m
*/ 
matrix GetParams(x,k,r,d,m){
	// create some variables that will get filled with transformed data
	Z0 = 0; Z1 = 0; Z2 = 0;
	TransformData(x,k,d,m, &Z0, &Z1, &Z2);
	T  = rows(Z0);
	p  = cols(Z0);
	p1 = cols(Z1);
	
	// if there are no lags then don't need to run FWL regression to 
	// concentrate out gammaHat
	if(k==0){
		R0 = Z0;
		R1 = Z1;
	}
	// lags are included, so use FWL
	else{
		Z2TZ2I = luinv(Z2'Z2);
		R0 = Z0 - Z2*Z2TZ2I*Z2'Z0;
		R1 = Z1 - Z2*Z2TZ2I*Z2'Z1;
	}

	S00 = R0'R0:/T;
	S01 = R0'R1:/T;
	S10 = R1'R0:/T;
	S11 = R1'R1:/T;
	if (r == 0){
		betaHat = J(p,1,0);
		betaStar = J(p1,1,0);
		alphaHat = J(p,1,0);
		PiHat = J(p,p,0);
		OmegaHat = S00;
  		if (m==0) rhoHat = J(p,1,.);  
		if (m==1){ 
		   rhoHat = J(p,1,0);
		}
	}
	if(0<r && r<p){
		// Reduced rank regression
		RRmatrix = luinv(S11)*S10*luinv(S00)*S01;
		Vtemp = .;
		D = .;
		// get the eigenvalues and eigen vectors
		eigensystem(RRmatrix,Vtemp,D);
		
		// sort the eigenvectors by their eigenvalues
		Vtemp = sort((Vtemp' , D'), p1+1);
		
		V = Vtemp[1..p1,.];
		
		/* betastar is a matrix of eigenvectors that correspond
		  to the r smallest eigenvalues
		  Note also, that if a constant is included, it is in the last
		  row of the beta start matrix
		*/
		betaStar = V[1..p1, p1..p1-r+1];
		
		// solve for the other parameters
		alphaHat = S01*betaStar*luinv(betaStar'*S11*betaStar);
		OmegaHat = S00 - alphaHat*betaStar'*S11*betaStar*alphaHat';	
		betaHat = betaStar[1..p, 1..r];
		PiHat = alphaHat*betaHat';

		G = luinv(betaHat[1..r,1..r]);
		betaHat = betaHat*G;
		betaStar = betaStar*G;
		alphaHat = alphaHat*luinv(G)';
		if (m==0) rhoHat = J(p,r,.);		
		if (m==1){ 
		   rhoHat = J(p,1,1)*betaStar[p1,.];
		}
 	}
	if(r == p){
		// if the system is full rank beta=pi and alpha = I_p
		V = S01*luinv(S11);
		betaHat = V[.,1..p]';
		alphaHat = I(p);
		PiHat = betaHat;
		OmegaHat = S00 - S01*luinv(S11)*S10;
		betaStar = betaHat;
		if (m==0) rhoHat = J(p,r,.);
		if (m==1){ 
			rhoHat = V[.,p1]';
			betaStar = betaHat \ rhoHat ;
			rhoHat = J(p,1,1) * rhoHat;
		}
	}
	
	PiStar = alphaHat*betaStar';
	
	if (k == 0)
		GammaHat = J(p,p,.);
	else{
		// REVERSE FWL to get GammaHat's
		GammaHat = Z2TZ2I*Z2'*(Z0 - Z1*PiStar');
		GammaHat = GammaHat';
	}
		
	parameters = OmegaHat , alphaHat , betaHat , rhoHat , GammaHat;
	return(parameters)
}


/* 
   this function calculates the likelihood
*/
matrix FCVARlike(x, k, r, d, m){
	T = rows(x);
	p = cols(x);
	estimates = GetParams(x, k, r, d, m);
	OmegaHat = estimates[1..p, 1..p];
	like = - (T*p)/2*( log(2*pi()) + 1)  - T/2*log(det(OmegaHat));
	return(like);
}

/* 
   This function calculates the likelihood over a grid of d's
   It is a useful first step in estimation because the estimation procedure
   can be sensitive to starting values when the likelihood function
   has multiple peaks (i.e. local max's), probelm explored by Federico Carlini (CREATES)   
*/
matrix LikeGrid(x,k,r,m,dMin,dMax,dStep){
	dGrid = range(dMin, dMax, dStep) 
	N     = rows(dGrid);
	likeD = J(N,1,0);
	
	for(i=1; i<=N; i++){
		d = dGrid[i,1];
		likeD[i,1] = Re( FCVARlike(x, k, r, d, m) );
	}
	
	iMax = .;
	maxMatrix = .;
	
	maxindex(likeD,1,iMax, maxMatrix);

	dStar = dGrid[iMax, 1];
	
	return(dStar);
}




/* 
   This is the way that function needs to be specified in Mata for the optimizer
    d is the paramater, (x,k,r,m) are arguments, ll is the returned value
    The mata optimizer does not allow for inequality constrains and it also
    doesn't allow restrictions on the bounds of the parameter space. To get around
    this, I impose the bounds with "if" statements.    
*/
function loglik(todo, d, x, k, r, m, dmin, dmax, ll, g, H){
	
	if(d<dmin){
		d = dmin;
	}
	if(d>dmax){
		d = dmax;
	}
	
	ll = Re( FCVARlike(x,k,r,d,m) )

}

// This function optimizes the above function given a starting value
scalar SolveModel(x, k, r, m, dstart, dmin, dmax){
	s=optimize_init()
	optimize_init_evaluator(s,&loglik())
	optimize_init_params(s,dstart)
	optimize_init_argument(s,1,x)
	optimize_init_argument(s,2,k)
	optimize_init_argument(s,3,r)
	optimize_init_argument(s,4,m)
	optimize_init_argument(s,5,dmin)
	optimize_init_argument(s,6,dmax)
	optimize_init_conv_maxiter(s, 500)
	dHat = optimize(s)
	return(dHat)
}

/* 
  this function calculates the residuals from the regression
*/
real matrix GetResiduals(x, k, r, d, m, alpha, beta, rho, Gamma){
	Z0 = .; 
	Z1 = .; 
	Z2 = .;
	TransformData(x,k,d,m, &Z0, &Z1, &Z2);
	epsilon = Z0;
	
	if (r > 0){
		if(m==1){
			epsilon = epsilon - Z1*( beta \ rho )* alpha';
		}
		if(m==0){
			epsilon = epsilon - Z1*( beta )* alpha';
		}
	}
	if (k > 0){
		epsilon = epsilon - Z2*Gamma';
	}
	
	epsilon = Re (epsilon) ;
	
	return(epsilon)
}

/* 
   this function returns the unconcentrated likelihood given parameters, it is
   used to obtain the Hessian matrix 
*/
real matrix FullFCVARlike(x, k, r, d, m, alpha, beta, rho, Gamma){
	T = rows(x);
	p = cols(x);
	epsilon = GetResiduals(x, k, r, d, m, alpha, beta, rho, Gamma);
	OmegaHat = epsilon'epsilon;
	like = - (T*p)/2*( log(2*pi()) + 1)  - T/2*log(det(OmegaHat));
	return(like);
}

/* 
   this function calculates the hessian matrix by numerically calculating 
   the second derivative for the parameters with known distributions
   Note: this includes alpha, gamma, and d
*/
real matrix FCVARhess(x, k, r, dHat, m, estimates){
	p = cols(x);
	
	nCols = cols(estimates)
	
	// set up an index to recover estimates
	i = r
	if(r==0){
		i = 1
	}

	// recover the estimates
	alphaHat = estimates[1..p, (p+1)..(p+1)+i-1];
	betaHat  = estimates[1..p, (p+1+i)..(p+1+i)+i-1];
	rhoHat   = estimates[1, (p+1+2*i)..(p+1+2*i)+i-1];
	gammaHat = estimates[1..p, (p+1+3*i)..nCols];
	
	// initialized perturbed variables
	alpha = .;
	gamma = .;
	d = .;
	
	/* create a vector of parameters for which the Hessian
	   matrix will be estimated
	*/
	pars = dHat;

	if(r>0){
		pars = pars \ vec( alphaHat );
	}
	if(k>0){
		pars = pars  \ vec( gammaHat );
	}
	
	nPars = length(pars);
	
	hessian = J(nPars,nPars,0);

	/* 
	   delta is the amount by which parameters are changed to calcualte
	   derivatives
	*/ 
	delta = 10^-4;
	
	/*
	  The following loop works in this way: like1, like2, like3, like4 
	  all measure the change in likelihood resulting from a delta shift
	  of the parameters. i.e. (like1 - like2) \ (2*delta) gives an approximation
	  of the first derivative. The second derivative is approximated by:
	  [(like1 - like2) \ (2*delta) - (like1 - like2) \ (2*delta)] \ (2*delta) 
	  Since he Hessian is symmetric, only half of the matrix needs to calculated
	*/
	
	for (s = 1; s<=nPars; s++){
	
	    for (j = 1; j<=s; j++){
	    
			pars1 = J(nPars,1,0);
			
			pars1[s] = delta;
			pars1[j] = pars1[j]  + delta;
			
			AdjustPars( pars1, &alpha, &d, &gamma, alphaHat, dHat, gammaHat, k, r, p);
			like1 = FullFCVARlike(x, k, r, d, m, alpha, betaHat, rhoHat, gamma);
			
			pars1[s] = -delta;
			pars1[j] = pars1[j]  + delta;
			AdjustPars( pars1, &alpha, &d, &gamma, alphaHat, dHat, gammaHat, k, r, p);
			like2 = FullFCVARlike(x, k, r, d, m, alpha, betaHat, rhoHat, gamma);			
			
			pars1[s] = delta;
			pars1[j] = pars1[j] -delta;
			AdjustPars( pars1, &alpha, &d, &gamma, alphaHat, dHat, gammaHat, k, r, p);
			like3 = FullFCVARlike(x, k, r, d, m, alpha, betaHat, rhoHat, gamma);			
			
			pars1[s] = -delta;
			pars1[j] = pars1[j] -delta;
			AdjustPars( pars1, &alpha, &d, &gamma, alphaHat, dHat, gammaHat, k, r, p);
			like4 = FullFCVARlike(x, k, r, d, m, alpha, betaHat, rhoHat, gamma);
			
			hessian[s,j] = (like1 - like2 - like3 + like4) / (4*delta^2);
			
			hessian[j,s] = hessian[s,j];
			
	    }
	    
	}
	
	return(hessian);
	
}

/* 
    This function uses pointers to adjust the alpha, d, and gamma parameters 
    adjusted by the delta shifts for the hessian matrix
*/
void AdjustPars( parameters, alpha, d, gamma, alphaHat, dHat, gammaHat, k, r, p){

	*d = dHat + parameters[1]
	N = length(parameters);
	if(r>0){
		*alpha = alphaHat + colshape( vec( parameters[2..p*r+1] ), p )';
	}
	if(k>0){
		*gamma = gammaHat + rowshape( vec( parameters[p*r+2..N] ), p*k )';
	}
}

/* 
   This function calls the LikeGrid function to get the starting value
   for a user specified grid of d's. It then estimates the model based on this
   approximate value as a starting point. 
   It then prints all of the coefficient estimates.
*/
scalar SolveModelGrid(x, k, r, m, dMin, dMax, dStep){
	// STEP 1: GRID SEARCH
	dStart = LikeGrid(x, k, r, m, dMin, dMax, dStep);
	
	// set lower bound
	dLB = dStart - 0.1;
	
	// set upper bound
	dUB = dStart + 0.1;
	

	// STEP 2: NUMERICAL REFINEMENT
	s=optimize_init()
	optimize_init_evaluator(s, &loglik())
	optimize_init_params(s, dStart)
	optimize_init_argument(s, 1, x)
	optimize_init_argument(s, 2, k)
	optimize_init_argument(s, 3, r)
	optimize_init_argument(s, 4, m)
	optimize_init_argument(s, 5, dLB);
	optimize_init_argument(s, 6, dUB);
	optimize_init_conv_maxiter(s, 30)
	
	// RECOVER THE RESULTS
	dHat = optimize(s)
	
	ll = optimize_result_value(s);
	
	ll = Re( ll );
	
	st_matrix("logl", ll);
	
	estimates = GetParams(x, k, r, dHat, m);
	
	estimates = Re( estimates );
	
	p = cols(x);
	
	T = rows(x);
	
	nCols = cols(estimates);
	
	
	// the rest of the function prints the estimation results
	i = r
	if(r==0){
		i = 1
	}
	OmegaHat = estimates[1..p, 1..p];
	alphaHat = estimates[1..p, (p+1)..(p+1)+i-1];
	betaHat  = estimates[1..p, (p+1+i)..(p+1+i)+i-1];
	rhoHat   = estimates[1, (p+1+2*i)..(p+1+2*i)+i-1];
	gammaHat = estimates[1..p, (p+1+3*i)..nCols];
	
	
	H = FCVARhess(x, k, r, dHat, m, estimates);
	SE = sqrt(-diagonal(luinv(H)));
	
	nPars = length(SE);
	

//	estimates = estimates[1..p, (p+1)..(p+1)+i-1], estimates[1..p, (p+1+3*i)..nCols];
//	estimates = dHat \ vec(estimates);
	
	printf("---------------------------------------------------------\n");
	printf("Estimation Results:\n");
	printf("---------------------------------------------------------\n");
	printf("Log-likelihood: %g\n", ll);
	printf("Number of variables: %g\n", p);
	printf("Number of observations: %g\n", T);
	printf("Number of lags: %g\n", k);
	printf("Rank: %g\n", r);
	printf("---------------------------------------------------------\n");
	printf("Fractional Differencing Parameter\n");
	printf("{txt}{space 13}{c |}      Coef.    Std. Err.\n");
	printf("{hline 13}{c +}{hline 24}\n");
        printf("{txt}%12s {c |} {res}%10.0g  %10.0g\n",
                                                "d", dHat, SE[1] )
	

	
	if(r>0){
	piHat = alphaHat*betaHat';

	alphaSE = colshape( vec( SE[2..p*r+1] ), p )';
	printf("{txt}Cointegrating Vector(s) \n");
	printf("{txt}{space 14}{c |} {res}%10.0g", 1);
		for(s = 2; s <=r; s++){
				printf("  %10.0g", s);
			}
			printf("\n");
	printf("{hline 14}{c +}{hline 12}");
		for(s = 2; s <=r; s++){
			printf("{hline 12}");
		}
		printf("\n");

		for(j = 1; j <=p; j++){
			printf("{txt}%12s%g {c |} {res}%10.0g", "beta",j, betaHat[j,1])
			for(s = 2; s <=r; s++){
				printf("  %10.0g", betaHat[j,s]);
			}
			printf("\n");
		}
	if(m==1){
		printf("{txt}%12s  {c |} {res}%10.0g", "rho", rhoHat[1])
		for(j = 2; j <=r; j++){
			printf("  %10.0g", rhoHat[j]);
			
		}
		printf("\n");
	
	}
	printf("{txt}Adjustment Matrix \n");
	printf("{txt}{space 14}{c |} {res}%10.0g", 1);
		for(s = 2; s <=r; s++){
			printf("  %10.0g", s);
		}
		printf("\n");
	printf("{hline 14}{c +}{hline 12}");
		for(s = 2; s <=r; s++){
			printf("{hline 12}");
		}
		printf("\n");
		for(j = 1; j <=p; j++){
			printf("{txt}%12s%g {c |} {res}%10.0g", "alpha",j, alphaHat[j,1])
			for(s = 2; s <=r; s++){
				printf("  %10.0g", alphaHat[j,s]);
			}
			printf("\n")
			printf("{txt}%12s  {c |} {res}%10.0g", "SE", alphaSE[j,1])
			for(s = 2; s <=r; s++){
				printf("  %10.0g", alphaSE[j,s]);
			}
			
			printf("\n");
		}
	printf("{txt}PI matrix\n")
			printf("{txt}{space 14}{c |} {res}%10.0g", 1);
			for(s = 2; s <=p; s++){
				printf("  %10.0g", s);
			}
			printf("\n");
			printf("{hline 14}{c +}{hline 12}");
			for(s = 2; s <=p; s++){
				printf("{hline 12}");
			}
			printf("\n");

		for(j = 1; j <=p; j++){
			printf("{txt}%12s%g {c |} {res}%10.0g", "Pi",j, piHat[j,1])
			for(s = 2; s <=p; s++){
				printf("  %10.0g", piHat[j,s]);
			}
			printf("\n");
		}
	}
	
	if(k>0){
	gammaSE = rowshape( vec( SE[p*r+2..nPars] ), p*k )';
		printf("{txt}Coefficients on Lags \n");
		for(l = 1; l <=k; l++){
			
			printf("{txt}Lag %g\n", l)
				printf("{txt}{space 14}{c |} {res}%10.0g", 1);
			for(s = 2; s <=p; s++){
				printf("  %10.0g", s);
			}
			printf("\n");
			printf("{hline 14}{c +}{hline 12}");
			for(s = 2; s <=p; s++){
				printf("{hline 12}");
			}
			printf("\n");

		for(j = 1; j <=p; j++){
			printf("{txt}%12s%g {c |} {res}%10.0g", "Gamma",j, gammaHat[j,p*(l-1)+1])
			for(s = 2; s <=p; s++){
				printf("  %10.0g", gammaHat[j,p*(l-1)+s]);
			}
			printf("\n")
			printf("{txt}%12s  {c |} {res}%10.0g", "SE", gammaSE[j,p*(l-1)+1])
			for(s = 2; s <=p; s++){
				printf("  %10.0g", gammaSE[j,p*(l-1)+s]);
			}
			
			printf("\n");
		}
		}
	}
		
	return(dHat)
}
end

