################################################################################
# 
# Library of Functions for FCVAR Estimation
# 
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 11, 2020
# 
################################################################################
# 
# Written by Michal Popiel and Morten Nielsen (This version 01.10.2020)
#
# DESCRIPTION: A library of functions called by the FCVAR estimation function. 
# 
################################################################################


################################################################################
# Load packages
################################################################################

# Load packages. 
# install.packages('name_of_package')
# library(name_of_package)



################################################################################
# Define functions.
################################################################################


################################################################################
# Define function to calculate cointegrating parameters from the data. 
################################################################################
# 
# function [ estimates ] <- GetParams(x, k, r, db, opt)
# Written by Michal Popiel and Morten Nielsen (This version 04.13.2016)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function uses FWL and reduced rank regression to obtain
# the estimates of Alpha, Beta, Rho, Pi, Gamma, and Omega
# 
# Input <- x   (matrix of variables to be included in the system)
#         k   (number of lags)
#         r   (number of cointegrating vectors)
#         db  (value of d and b)
#         opt (object containing the estimation options)
# Output <- estimates (Matlab structure containing the following)
#          - estimates.db (taken directly from the input)
#          - estimates.alphaHat
#          - estimates.betaHat
#          - estimates.rhoHat
#          - estimates.piHat
#          - estimates.OmegaHat
#          - estimates.GammaHat ( p x kp matrix [GammaHat1,...,GammaHatk])
#
################################################################################


GetParams <- function(x, k, r, db, opt) {
  
  # print('Made it to GetParams!')
  # print('k = ')
  # print(k)
  # print('r = ')
  # print(r)
  # print('db = ')
  # print(db)
  
  Z_array <- TransformData(x, k, db, opt)
  Z0 <- Z_array$Z0
  Z1 <- Z_array$Z1
  Z2 <- Z_array$Z2
  Z3 <- Z_array$Z3
  
  T <- nrow(Z0)
  p <- ncol(Z0)
  p1 <- ncol(Z1)
  
  
  #--------------------------------------------------------------------------------
  # Concentrate out the unrestricted constant if present
  #--------------------------------------------------------------------------------
  
  if(opt$unrConstant) {
    # Note Z3*inv(Z3*Z3)*Z3 is just a matrix of ones / T
    Z0hat <- Z0 - Z3 %*% ( solve(t(Z3) %*% Z3) %*% t(Z3) %*% Z0 )
    Z1hat <- Z1 - Z3 %*% ( solve(t(Z3) %*% Z3) %*% t(Z3) %*% Z1 )
    if(k > 0) {
      Z2hat <- Z2 - Z3 %*% ( solve(t(Z3) %*% Z3) %*% t(Z3) %*% Z2 )
    } else {
      Z2hat <- Z2
    }
    
  } else {
    Z0hat <- Z0
    Z1hat <- Z1
    Z2hat <- Z2
  }
  
  #--------------------------------------------------------------------------------
  # FWL Regressions
  #--------------------------------------------------------------------------------
  
  if (k == 0) {
    # No lags, so there are no effects of Z2.
    R0 <- Z0hat
    R1 <- Z1hat
  } else {
    # Lags included: Obtain the residuals from FWL regressions.
    R0 <- Z0hat - Z2hat %*% ( solve(t(Z2hat) %*% Z2hat) %*% t(Z2hat) %*% Z0hat )
    R1 <- Z1hat - Z2hat %*% ( solve(t(Z2hat) %*% Z2hat) %*% t(Z2hat) %*% Z1hat )
  }
  
  
  # Calculate Sij matrices for reduced rank regression.
  S00 <- t(R0) %*% R0/T
  S01 <- t(R0) %*% R1/T
  S10 <- t(R1) %*% R0/T
  S11 <- t(R1) %*% R1/T
  
  # Calculate reduced rank estimate of Pi.
  if (r == 0) {
    
    betaHat <- NULL
    betaStar <- NULL
    alphaHat <- NULL
    PiHat <- NULL
    rhoHat <- NULL
    OmegaHat <- S00
    
  } else if (( r > 0 ) & ( r < p )) {
    
    # For d = startVal, need to match this:
    # >> [V,D] = eig( inv(S11)*S10*inv(S00)*S01 );
    # >> V
    # 
    # V =
    #   
    #   -0.9441   -0.8883    0.3380
    # -0.1072    0.0466    0.0062
    # -0.3117    0.4569    0.9411
    # 
    # >> D
    # 
    # D =
    #   
    #   0.0506         0         0
    # 0    0.0257         0
    # 0         0    0.0015
    # 
    # >>
    #   
    
    
    # Clean this up:
    eig_out <- eigen( solve(S11) %*% S10 %*% solve(S00) %*% S01 )
    D <- eig_out$values # Only the vector, not the diagonal matrix. 
    V1 <- eig_out$vectors
    
    # V <- sortrows( [ V diag(D) ], p1+1 )
    V2 <- t( cbind(t(V1), D)[order(D), ] )
    
    V <- V2[1:p1, ]
    # betaStar <- V( 1:p1, p1 : -1 : p1-r+1 )
    betaStar <- matrix(V[ 1:p1, seq(p1, p1-r+1, by <- -1) ], 
                       nrow = p1, ncol = p1-r-1)
    
    
    
    
    # If either alpha or beta is restricted, then the likelihood is
    #   maximized subject to the constraints. This section of the code
    #   replaces the call to the switching algorithm in the previous
    #   version.
    if(!is.null(opt$R_Alpha) | !is.null(opt$R_Beta) ) {
      
      
      # [ betaStar, alphaHat, OmegaHat ]...
      # <- RstrctOptm_Switch(betaStar, S00, S01, S11, T, p, opt)
      switched_mats <- RstrctOptm_Switch(betaStar, S00, S01, S11, T, p, opt)
      betaStar <- switched_mats$betaStar
      alphaHat <- switched_mats$alphaHat
      OmegaHat <- switched_mats$OmegaHat
      
      betaHat <- betaStar[1:p, 1:r]
      PiHat <- alphaHat %*% t(betaHat) 
      
      
      
    } else {
      
      
      # Otherwise, alpha and beta are unrestricted, but unidentified.
      alphaHat <- S01 %*% betaStar %*% solve(t(betaStar) %*% S11 %*% betaStar) 
      OmegaHat <- S00 - alphaHat %*% t(betaStar) %*% S11 %*% betaStar %*% t(alphaHat)
      betaHat <- matrix(betaStar[1:p, 1:r], nrow = p, ncol = r)
      PiHat <- alphaHat %*% t(betaHat) 
      # Transform betaHat and alphaHat to identify beta.
      #   The G matrix is used to identify beta by filling the first
      #   rxr block of betaHat with an identity matrix.
      G <- solve(betaHat[1:r,1:r])
      betaHat <- betaHat %*% G
      betaStar <- betaStar %*% G
      # alphaHat is post multiplied by G^-1 so that Pi= a(G^-1)Gb <- ab
      alphaHat <- alphaHat %*% t(solve(G)) 
      
    }
    
    
    
    # Extract coefficient vector for restricted constant model if required.
    if (opt$rConstant) {
      rhoHat <- betaStar[p1, ]
    }
    else {
      rhoHat <- NULL
    }
    
    
  } else {
    
    # (r <- p) and do no need reduced rank regression
    V <- S01 %*% solve(S11)
    betaHat <- t(V[,1:p]) 
    # For restriced constant, rho <- last column of V.
    alphaHat <- diag(p)
    PiHat <- betaHat
    OmegaHat <- S00 - S01 %*% solve(S11) %*% S10
    
    # Extract coefficient vector for restricted constant model if required.
    if (opt$rConstant) {
      rhoHat <- t(V[ , p1]) 
    } else {
      rhoHat <- NULL
    }
    
    
    # Check conformability: 
    betaStar <- cbind(betaHat, rhoHat) 
    
  }
  
  
  
  # Calculate PiStar independently of how betaStar was estimated.
  PiStar <- alphaHat %*% t(betaStar) 
  
  # Perform OLS regressions to calculate unrestricted constant and
  #   GammaHat matrices if necessary.
  xiHat <- NULL
  if (k == 0) {
    GammaHat <- NULL
  } else {
    
    if(r > 0) {
      GammaHat <- t( solve(t(Z2hat) %*% Z2hat) %*% t(Z2hat) %*% (Z0hat - Z1hat %*% t(PiStar)) )
    } else {
      GammaHat <- t( solve(t(Z2hat) %*% Z2hat) %*% t(Z2hat) %*% Z0hat ) 
    }
    
  }
  
  
  
  if(opt$unrConstant==1) {
    if(k>0) {
      if(r>0) {
        xiHat <- t( solve(t(Z3) %*% Z3) %*% t(Z3) %*% (Z0 - Z1 %*% t(PiStar) - Z2 %*% t(GammaHat)) )  
      }
      else {
        xiHat <- t( solve(t(Z3) %*% Z3) %*% t(Z3) %*% (Z0 - Z2 %*% t(GammaHat)) ) 
      }
      
    } else {
      if(r>0) {
        xiHat <- t( solve(t(Z3) %*% Z3) %*% t(Z3) %*% (Z0 - Z1 %*% t(PiStar)) ) 
      }
      else {
        xiHat <- t( solve(t(Z3) %*% Z3) %*% t(Z3) %*% Z0 )  
      }
      
    }
    
  }
  
  
  
  
  #--------------------------------------------------------------------------------
  # Return the results in a list
  #--------------------------------------------------------------------------------
  
  estimates <- list(
    
    db = db, 
    alphaHat = alphaHat, 
    betaHat = betaHat, 
    rhoHat = rhoHat, 
    PiHat = PiHat, 
    OmegaHat = OmegaHat, 
    GammaHat = GammaHat, 
    xiHat = xiHat
    
  )
  
  
  return(estimates)
  
}






################################################################################
# Define function to calculate the likelihood for the unconstrained model. 
################################################################################
# 
# function [ like ] <- FCVARlikeMu(y, db, mu, k, r, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# 
# DESCRIPTION: This function evaluates the likelihood for a given set of
# 	parameter values. It is used by the LikeGrid() function to numerically
# 	optimize over the level parameter for given values of the fractional
# 	parameters.
# 
# Input <- y  (matrix of variables to be included in the system)
#         db (fractional parameters d,b)
#         mu (level parameter)
#         k  (number of lags)
#         r  (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- like (log-likelihood evaluated at specified parameter values)
# 
################################################################################



FCVARlikeMu <- function(mu, y, db, k, r, opt) {
  
  
  t <- nrow(y)
  x <- y - matrix(1, nrow = t, ncol = 1) %*% mu
  
  # Obtain concentrated parameter estimates.
  estimates <- GetParams(x, k, r, db, opt)
  
  # Calculate value of likelihood function.
  T <- t - opt$N
  p <- ncol(x)
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates$OmegaHat))
  
  
  
  return(like)
}




################################################################################
# Define function to calculate the likelihood for the constrained model. 
################################################################################
# 
# function [ like ] <- FCVARlike(x, params, k, r, opt)
# Written by Michal Popiel and Morten Nielsen (This version 11.10.2014)
# 
# DESCRIPTION: This function adjusts the variables with the level parameter,
# 	if present, and returns the log-likelihood given d,b.
# 
# Input <- x   (matrix of variables to be included in the system)
#         params (a vector of parameters d,b, and mu (if option selected))
#         k   (number of lags)
#         r   (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- like (concentrated log-likelihood evaluated at given parameters)
# 
################################################################################



FCVARlike <- function(params, x, k, r, opt) {
  
  # global estimatesTEMP
  
  T <- nrow(x)
  p <- ncol(x)
  
  # print('params = ')
  # print(params)
  
  
  # If there is one linear restriction imposed, optimization is over phi,
  #  so translate from phi to (d,b), otherwise, take parameters as
  #  given.
  if(nrow(opt$R_psi) == 1) {
    H <- null(opt$R_psi)
    h <- t(opt$R_psi) %*% solve(opt$R_psi %*% t(opt$R_psi)) %*% opt$r_psi       
    db <- H*params[1] + h
    d <- db[1]
    b <- db[2]
    if (opt$levelParam) {
      mu <- params[2:length(params)]
    }
    
  } else {
    d <- params[1]
    b <- params[2]    
    if (opt$levelParam) {
      mu <- params[3:length(params)]
    }
    
  }
  
  
  # Modify the data by a mu level.
  if (opt$levelParam) {
    y <- x - matrix(1, nrow = T, ncol = p) %*% diag(mu)
  } else {
    y <- x
  }
  
  
  # If d=b is not imposed via opt$restrictDB, but d>=b is required in
  #	opt$constrained, then impose the inequality here.
  if(opt$constrained) {
    b <- min(b, d)
  }
  
  
  # Obtain concentrated parameter estimates$
  estimates <- GetParams(y, k, r, c(d, b), opt)
  
  # Storing the estimates in a global structure here allows us to skip a
  #   call to GetParams after optimization to recover the coefficient
  #   estimates
  estimatesTEMP <- estimates
  # If level parameter is present, they are the last p parameters in the
  #   params vector
  if (opt$levelParam) {
    estimatesTEMP$muHat <- mu
  } else {
    estimatesTEMP$muHat <- NULL
  }
  
  
  # Calculate value of likelihood function.
  T <- nrow(y) - opt$N
  p <- ncol(y)
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates$OmegaHat))
  
  # Assign the value for the global variable estimatesTEMP.
  # Might adjust this later but it will work for now. 
  # print('Search list in FCVARlike():')
  # print(search())
  assign("estimatesTEMP", estimatesTEMP, envir = .GlobalEnv)
  
  return(like)
}




################################################################################
# Define function to calculate the likelihood, 
# evaluated at the parameters provided as inputs. 
################################################################################
# 
# function [ like ] <- FullFCVARlike(x, k, r, coeffs, beta, rho, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function returns the value of the log-likelihood
# 	evaluated at the parameters provided as inputs.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         coeffs (Matlab structure of coefficients)
#         beta (value of beta)
#         rho  (value of rho)
#         opt  (object containing the estimation options)
# Output <- like (value of the log likelihood)
# 
################################################################################

FullFCVARlike <- function(x, k, r, coeffs, beta, rho, opt) {
  
  
  # Add betaHat and rhoHat to the coefficients to get residuals because they
  #   are not used in the Hessian calculation and are missing from the
  #   structure coeffs
  coeffs$betaHat <- beta
  coeffs$rhoHat <- rho
  
  # Obtain residuals.
  epsilon <- GetResiduals(x, k, r, coeffs, opt)
  
  
  # Calculate value of likelihood function.
  T <- nrow(x) - opt$N
  p <- ncol(x)
  OmegaHat <- t(epsilon) %*% epsilon/T
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(OmegaHat)) 
  
  return(like)
}



################################################################################
# Define function to transform the raw data by fractional differencing.
################################################################################
# 
# function [ Z0, Z1, Z2, Z3 ] <- TransformData(x, k, db, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (May 24, 2013)
# 
# DESCRIPTION: Returns the transformed data required for regression and
# 	reduced rank regression.
# 
# Input <- x   (matrix of variables to be included in the system)
#         k   (number of lags)
#         db  (fractional differencing parameters d and b)
#         opt (object containing the estimation options)
# Output <- Z0, Z1, Z2, and Z3 of transformed data.
# 
# Calls the function FracDiff(x, d) and Lbk(x, b, k).
# 
################################################################################

TransformData <- function(x, k, db, opt) {
  
  
  # Number of initial values and sample size.
  N <- opt$N
  T <- nrow(x) - N
  
  # Extract parameters from input.
  d <- db[1]
  b <- db[2]
  
  # Transform data as required.
  Z0 <- FracDiff(x, d)
  
  Z1 <- x
  # Add a column with ones if model includes a restricted constant term.
  if(opt$rConstant) {
    Z1 <- cbind(x, matrix(1, nrow = N+T, ncol = 1))
  }
  
  Z1 <- FracDiff(  Lbk( Z1 , b, 1)  ,  d - b )
  
  Z2 <- FracDiff(  Lbk( x , b, k)  , d)
  
  # Z3 contains the unrestricted deterministics
  Z3 <- NULL
  # # Add a column with ones if model includes a unrestricted constant term.
  if(opt$unrConstant) {
    Z3 <- matrix(1, nrow = T, ncol = 1)
  }
  
  
  # Trim off initial values.
  Z0 <- Z0[(N+1):nrow(Z0), ]
  Z1 <- Z1[(N+1):nrow(Z1), ]
  if(k > 0) {
    Z2 <- Z2[(N+1):nrow(Z2), ]
  }
  
  
  # Return all arrays together as a list. 
  Z_array <- list(
    Z0 = Z0,
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3
  )
  
  return(Z_array)
}


################################################################################
# Define function to calculate residuals from given parameter values.
################################################################################
# 
# function [ epsilon ] <- GetResiduals(x, k, r, coeffs, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function calculates the model residuals.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         coeffs (Matlab structure of coefficients)
#         opt (object containing the estimation options)
# Output <- epsilon (matrix of residuals from model estimation evaluated at
#                     the parameter estimates specified in coeffs) 
# 
################################################################################

GetResiduals <- function(x, k, r, coeffs, opt) {
  
  
  # If level parameter is included, the data must be shifted before
  #   calculating the residuals:
  if (opt$levelParam) {
    T <- size(x,1)
    y <- x - matrix(1, nrow = T, ncol = 1) %*% coeffs$muHat
  } else {
    y <- x
  }
  
  #--------------------------------------------------------------------------------
  # Transform data 
  #--------------------------------------------------------------------------------
  # [ Z0, Z1, Z2, Z3 ] <- TransformData(y, k, coeffs$db, opt)
  Z_array <- TransformData(y, k, coeffs$db, opt)
  Z0 <- Z_array$Z0
  Z1 <- Z_array$Z1
  Z2 <- Z_array$Z2
  Z3 <- Z_array$Z3
  
  #--------------------------------------------------------------------------------
  # Calculate residuals 
  #--------------------------------------------------------------------------------
  epsilon <- Z0
  
  if (r > 0) {
    epsilon <- epsilon - 
      Z1 %*% cbind(coeffs$betaHat, coeffs$rhoHat) %*% t(coeffs$alphaHat)
  }
  
  
  if (k > 0) {
  epsilon <- epsilon - Z2 %*% t(coeffs$GammaHat)
  }
  
  
  if (opt$unrConstant) {
    epsilon <- epsilon - Z3 %*% t(coeffs$xiHat)
  }
  
  
  
  return(epsilon)
}
end




################################################################################
# Define function to calculate lag polynomial in the fractional lag operator.
################################################################################
# 
# function [ Lbkx ] <- Lbk(x, b, k)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (May 24, 2013)
# 
# DESCRIPTION: Lbk(x, b, k) is a lag polynomial in the fractional lag operator.
# 
# Input <- x (vector or matrix of data)
#         b (scalar value at which to calculate the fractional lag)
#         k (number of lags)
# 
# Output <- matrix  [ Lb^1 x, Lb^2 x, ..., Lb^k x] where Lb <- 1 - (1-L)^b. The output 
#		matrix has the same number of rows as x but k times as many columns.
# 
# Calls the function FracDiff(x, d) 
# 
################################################################################


Lbk <- function(x, b, k) {
  
  
  p <- ncol(x)
  
  # Initialize output matrix.
  Lbkx <- NULL
  
  # For i <- 1, set first column of Lbkx <- Lb^1 x.
  if (k > 0) {
    bx <- FracDiff(x, b)
    Lbkx <- x - bx
  }
  
  
  if (k > 1) {
    for (i in 2:k ) {
      Lbkx <- cbind(Lbkx,  
                    ( Lbkx[ , p*(i-2)+1 : ncol(Lbkx)] - 
                        FracDiff(Lbkx[ , p*(i-2)+1 : ncol(Lbkx)], b) ))
      
    }  
    
  }
  
  
  
  return(Lbkx)
}



################################################################################
# Define function to fractionally difference a series. 
################################################################################
# 
# function [dx] = FracDiff(x,d)
# Andreas Noack Jensen & Morten Nielsen
# May 24, 2013
# 
# FracDiff(x,d) is a fractional differencing procedure based on the
# 	fast fractional difference algorithm of Jensen & Nielsen (2014, JTSA).
# 
# input = x (vector or matrix of data)
#         d (scalar value at which to calculate the fractional difference)
# 
# output = vector or matrix (1-L)^d x of same dimensions as x.
# 
################################################################################

FracDiff <- function(x, d) {
  
  
  T <- nrow(x)
  p <- ncol(x)
  
  # #--------------------------------------------------------------------------------
  # # Upgrade to the FFT version in a later version.
  # #--------------------------------------------------------------------------------
  # 
  # k <- matrix(seq(1, T-1), nrow = T-1, ncol = 1)
  # 
  # # NEXTPOW2(N) returns the first P such that 2.^P >= abs(N).  It is
  # #     often useful for finding the nearest power of two sequence
  # #     length for FFT operations.
  # NFFT <- 2^nextpow2(2*T-1)
  # 
  # # Array operation of the index of the series without the last element minus
  # # the order of integration+1, divided by that same series
  # b <- (k-d-1)/k
  # 
  # # cumulative product of that series modified in previous line
  # b <- c(1, cumprod(b))
  # 
  # # IFFT(X) is the inverse discrete Fourier transform of X.
  # #     IFFT(..., 'symmetric') causes IFFT to treat X as conjugate symmetric
  # #     along the active dimension.  This option is useful when X is not exactly
  # #     conjugate symmetric merely because of round-off error. 
  # # REPMAT Replicate and tile an array.
  # #     B = repmat(A,M,N) creates a large matrix B consisting of an M-by-N
  # #     tiling of copies of A. The size of B is [size(A,1)*M, size(A,2)*N].
  # #     The statement repmat(A,N) creates an N-by-N tiling.
  # dx <- ifft(repmat(fft(b, NFFT), 1, p) * fft(x, NFFT), 'symmetric')
  # 
  # dx <- dx[1:T, ]
  # 
  # 
  # #--------------------------------------------------------------------------------
  
  # For now, just use a wrapper for the diffseries() function in the fracdiff package. 
  # except that only differences the first column. 
  # This is based on the FFT version by Jensen and Nielsen (2014). 
  # x_diff <- diffseries(x, 0.8)
  # dx_MM <- data.frame(matrix(NA, nrow = nrow(x),
  #                         ncol = ncol(x)))
  # for (col_num in 1:length(colnames(x))) {
  #   dx_MM[, col_num] <- diffseries(x[, col_num], d)
  # }
  # head(dx_MM)
  # Problem: Their implementation does not match Morten's.
  
  # For now, I will go with the tried-and-tested "slow version".
  k <- matrix(seq(1, T-1), nrow = T-1, ncol = 1)
  b <- (k-d-1)/k
  b <- c(1, cumprod(b))

  xlead <- data.frame(matrix(0,
                             nrow = T,
                             ncol = ncol(x)))
  colnames(xlead) <- colnames(x)

  dx <- filter(rbind(xlead, x),
               filter = b,
               sides = 1)
  # Trim off the leading rows.
  dx <- dx[seq(T+1, 2*T), ]
  
  
  return(dx)
}



################################################################################
# Define function to counts the number of free parameters
################################################################################
# 
# function [ fp ] <- FreeParams(k, r, p, opt, rankJ)
# Written by Michal Popiel and Morten Nielsen (This version 05.22.2015)
# 
# DESCRIPTION: This function counts the number of free parameters based on
# 	the number of coefficients to estimate minus the total number of
# 	restrictions. When both alpha and beta are restricted, the rank condition
# 	is used to count the free parameters in those two variables.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- fp (number of free parameters)
# 
################################################################################


FreeParams <- function(k, r, p, opt, rankJ) {
  
  # ---- First count the number of parameters -------- %
  fDB <- 2 # updated by rDB below for the model d=b (1 fractional parameter)
  fpA <- p*r # alpha
  fpB <- p*r # beta
  fpG <- p*p*k # Gamma
  fpM <- p*opt$levelParam # mu
  fpRrh <- r*opt$rConstant # restricted constant
  fpUrh <- p*opt$unrConstant # unrestricted constant
  
  numParams <- fDB + fpA + fpB + fpG + fpM + fpRrh + fpUrh
  
  # ---- Next count the number of restrictions -------- %
  rDB <- nrow(opt$R_psi)
  
  
  if(isempty(opt$R_Beta) & isempty(opt$R_Alpha)) {
    # If Alpha or Beta are unrestricted, only an identification restriction
    #   is imposed
    rB <- r*r # identification restrictions (eye(r))
    numRest <- rDB + rB
    # ------ Free parameters is the difference between the two ---- %
    fp <- numParams - numRest
  }
  else {
    # If there are restrictions imposed on beta or alpha, we can use the rank
    #   condition calculated in the main estimation function to give us the
    #   number of free parameters in alpha, beta, and the restricted constant
    #   if it is being estimated:
    fpAlphaBetaRhoR <- rankJ
    fp <- fpAlphaBetaRhoR + (fDB + fpG + fpM + fpUrh) - rDB 
  }   
  
  
  return(fp)
}


################################################################################
# Define function to calculate the Hessian matrix. 
################################################################################
# 
# function [ hessian ] <- FCVARhess(x, k, r, coeffs, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# 
# DESCRIPTION: This function calculates the Hessian matrix of the
# 	log-likelihood numerically.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         coeffs (coefficient estimates around which estimation takes place)
#         opt (object containing the estimation options)
# Output <- hessian (matrix of second derivatives)
# 
################################################################################



FCVARhess <- function(x, k, r, coeffs, opt) {
  
  
  # Set dimensions of matrices.
  p <- ncol(x)
  
  # rhoHat and betaHat are not used in the Hessian calculation.
  rho <- coeffs$rhoHat
  beta <- coeffs$betaHat
  
  # Specify delta (increment for numerical derivatives).
  delta <- 10^(-4)
  # delta <- 10^(-8)
  
  # Convert the parameters to vector form.
  phi0 <- SEmat2vecU(coeffs, k, r, p, opt)
  
  # Calculate vector of increments.
  deltaPhi <- delta*matrix(1, nrow = nrow(phi0), ncol = ncol(phi0))
  
  nPhi <- length(phi0)
  
  # Initialize the Hessian matrix.
  hessian <- matrix(0, nrow = nPhi, ncol = nPhi)
  
  # The loops below evaluate the likelihood at second order incremental
  #   shifts in the parameters. 
  
  
  for (i in 1:nPhi) {
    
    
    for (j in 1:i) {
      
      
      # positive shift in both parameters.
      phi1 <- phi0 + deltaPhi*( (1:nPhi) == i ) + deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi1, k, r, p, opt )
      # calculate likelihood
      like1 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt )
      
      # negative shift in first parameter, positive shift in second
      # parameter. If same parameter, no shift.
      phi2 <- phi0 - deltaPhi*( (1:nPhi) == i ) + deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi2, k, r, p, opt )
      # calculate likelihood
      like2 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)
      
      # positive shift in first parameter, negative shift in second
      # parameter. If same parameter, no shift.
      phi3 <- phi0 + deltaPhi*( (1:nPhi) == i ) - deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi3, k, r, p, opt )
      # calculate likelihood
      like3 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)
      
      # negative shift in both parameters.
      phi4 <- phi0 - deltaPhi*( (1:nPhi) == i ) - deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi4, k, r, p, opt )
      # calculate likelihood
      like4 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)
      
      # The numerical approximation to the second derivative.
      hessian[i,j] <- ( like1 - like2 - like3 + like4 )/4/deltaPhi[i]/deltaPhi[j]
      
      # Hessian is symmetric.
      hessian[j,i] <- hessian[i,j]
      
      
    }
    
  }
  
  return(hessian)
}


################################################################################
# Define function to collect parameters into a vector
################################################################################
# 
# function [ param ] <- SEmat2vecU( coeffs, k, r, p , opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function transforms the model parameters in matrix
# 	form into a vector.
# 
# Input <- coeffs (Matlab structure of coefficients in their usual matrix form)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         p (number of variables in the system)
#         opt (object containing the estimation options)
# Output <- param (vector of parameters)
# 
################################################################################


SEmat2vecU <- function( coeffs, k, r, p , opt) {
  
  # If restriction d=b is imposed, only adjust d.
  if (opt$restrictDB) {
    param <- coeffs$db[1]
  }
  else {
    param <- coeffs$db
  }
  
  
  # Level parameter MuHat.
  if (opt$levelParam) {
    param <- cbind(param, matrix(coeffs$muHat, nrow = 1, ncol = p))
  }
  
  
  # Unrestricted constant.
  if (opt$unrConstant) {
    param <- cbind(param, matrix(coeffs$xiHat, nrow = 1, ncol = p))
  }
  
  
  # alphaHat
  if (r > 0) {
    param <- cbind(param, matrix( coeffs$alphaHat, nrow = 1, ncol = p*r ))
  }
  
  
  # GammaHat
  if (k > 0) {
    param <- cbind(param, matrix( coeffs$GammaHat, nrow = 1, ncol = p*p*k ))
  }
  
  
  return(param)
}



################################################################################
# Define function to extract parameters from a vector. 
################################################################################
# 
# function [ coeffs ] <- SEvec2matU( param, k, r, p, opt )
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function transforms the vectorized model parameters
# 	into matrices.
# 
# Input <- param (vector of parameters)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         p (number of variables in the system)
#         opt (object containing the estimation options)
# Output <- coeffs (Matlab structure of coefficients in their usual matrix form)
# 
################################################################################

SEvec2matU <- function( param, k, r, p, opt ) {
  
  
  if (opt$restrictDB) {
    coeffs$db <- matrix(c(param[1], param[1]), nrow = 1, ncol = 2)  # store d,b
    param <- param[2:length(param)]				# drop d,b from param
  } else {
    coeffs$db <- param[1:2]
    param <- param[3:length(param)]
  }
  
  
  if (opt$levelParam) {
    coeffs$muHat <- param[1:p]
    param <- param[(p+1):length(param)]
  }
  
  
  if (opt$unrConstant) {
    coeffs$xiHat <- matrix(param[1:p], nrow = p, ncol = 1)
    param <- param[(p+1):length(param)]
  }
  
  
  if (r > 0) {
    coeffs$alphaHat <- matrix( param[1:p*r], nrow = p, ncol = r)
    param <- param[(p*r+1):length(param)]
  } else {
    coeffs$alphaHat <- NULL
  }
  
  
  if (k > 0) {
    coeffs$GammaHat <- matrix( param[1 : length(param)], nrow = p, ncol = p*k)
  } else {
    coeffs$GammaHat <- NULL
  }
  
  
  return(coeffs)
}




################################################################################
# Define function to calculate the roots of the characteristic polynomial
################################################################################
# 
# function cPolyRoots <- CharPolyRoots(coeffs, opt, k, r, p)
# Written by Michal Popiel and Morten Nielsen (This version 12.07.2015)
# Based on Lee Morin & Morten Nielsen (May 31, 2013)
# 
# DESCRIPTION: CharPolyRoots calculates the roots of the 
#     characteristic polynomial and plots them with the unit circle 
#     transformed for the fractional model, see Johansen (2008).
# 
# input <- coeffs (Matlab structure of coefficients
#         opt (object containing the estimation options)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         p (number of variables in the system)
# 
# output <- complex vector cPolyRoots with the roots of the characteristic polynomial.
# 
# No dependencies.
# 
# Note: The roots are calculated from the companion form of the VAR, 
#       where the roots are given as the inverse eigenvalues of the 
#       coefficient matrix.
# 
################################################################################

CharPolyRoots <- function(coeffs, opt, k, r, p) {
  
  
  b <- coeffs$db[2]
  
  # First construct the coefficient matrix for the companion form of the VAR. 
  PiStar <- diag(p)
  if (r > 0) {
    PiStar <- PiStar + coeffs$alphaHat %*% t(coeffs$betaHat)
  }
  
  
  if (k > 0) {
    Gamma1 <- coeffs$GammaHat[ , 1 : p]
    PiStar <- PiStar + Gamma1
    for (i in 2:k) {
      
      Gammai <- coeffs$GammaHat[ , seq(((i-1)*p + 1), i*p)]
      GammaiMinus1 <- coeffs$GammaHat[ , seq(((i-2)*p + 1), (i-1)*p)]
      
      PiStar <- cbind(PiStar, (Gammai - GammaiMinus1))
      
    }
    
    Gammak <- coeffs$GammaHat[ , seq(((k-1)*p + 1), k*p)]
    PiStar <- cbind(PiStar, ( - Gammak ))
  }
  
  
  # Pad with an identity for the transition of the lagged variables.
  PiStar <- rbind(PiStar, 
                  cbind(diag(p*k), 
                        matrix(0, nrow = p*k, ncol = p*(k>0) )))
  
   
  # The roots are then the inverse eigenvalues of the matrix PiStar.
  cPolyRoots <- 1 / eigen(PiStar)$values
  cPolyRoots <- cPolyRoots[order(-Mod(cPolyRoots))]
  
  # Generate graph depending on the indicator plotRoots.
  if (opt$plotRoots) {
    # Now calculate the line for the transformed unit circle.
    # First do the negative half.
    unitCircle <- seq( pi, 0, by = - 0.001)
    psi <- - (pi - unitCircle)/2
    unitCircleX <- cos( - unitCircle)
    unitCircleY <- sin( - unitCircle)
    transformedUnitCircleX <- (1 - (2*cos(psi))^b*cos(b*psi))
    transformedUnitCircleY <- (    (2*cos(psi))^b*sin(b*psi))
    # Then do the positive half.
    unitCircle <- seq(0, pi, by = 0.001)
    psi <- (pi - unitCircle)/2
    unitCircleX <- c(unitCircleX, cos(unitCircle))
    unitCircleY <- c(unitCircleY, sin(unitCircle))
    transformedUnitCircleX <- c(transformedUnitCircleX, 1, 
                                (1 - (2*cos(psi))^b*cos(b*psi)))
    transformedUnitCircleY <- c(transformedUnitCircleY, 0, 
                                (    (2*cos(psi))^b*sin(b*psi)))
    
    # Plot the unit circle and its image under the mapping
    # along with the roots of the characterisitc polynomial.
    
    # Determine axes based on largest roots. 
    maxXYaxis <- max( c(transformedUnitCircleX, unitCircleX,
                        transformedUnitCircleY, unitCircleY) )
    minXYaxis <- min( c(transformedUnitCircleX, unitCircleX,
                        transformedUnitCircleY, unitCircleY) )
    maxXYaxis <- max( maxXYaxis, -minXYaxis )
    
    plot(transformedUnitCircleX, 
         transformedUnitCircleY, 
         main = c('Roots of the characteristic polynomial', 
                  'with the image of the unit circle'), 
         xlab = 'Real Part of Root', 
         ylab = 'Imaginary Part of Root', 
         xlim = 2*c(-maxXYaxis, maxXYaxis),
         ylim = 2*c(-maxXYaxis, maxXYaxis),
         type = 'l', 
         lwd = 3, 
         col = 'red')
    lines(unitCircleX, unitCircleY, lwd = 3, col = 'black')
    points(Re(cPolyRoots), Im(cPolyRoots), 
           pch = 16, col = 'blue')
    
    
    
  }
  
  
  return(cPolyRoots)
}






# Functions to make: 


#--------------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------------



################################################################################
# Define function to 
################################################################################
# 

# 
################################################################################




################################################################################
# End
################################################################################

