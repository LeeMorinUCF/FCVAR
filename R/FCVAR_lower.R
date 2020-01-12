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
# function [ estimates ] = GetParams(x, k, r, db, opt)
# Written by Michal Popiel and Morten Nielsen (This version 04.13.2016)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function uses FWL and reduced rank regression to obtain
# the estimates of Alpha, Beta, Rho, Pi, Gamma, and Omega
# 
# Input = x   (matrix of variables to be included in the system)
#         k   (number of lags)
#         r   (number of cointegrating vectors)
#         db  (value of d and b)
#         opt (object containing the estimation options)
# Output = estimates (Matlab structure containing the following)
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
  
  
  Z_array <- TransformData(x, k, db, opt)
  Z0 <- Z_array$mat_0
  Z1 <- Z_array$mat_1
  Z2 <- Z_array$mat_2
  Z3 <- Z_array$mat_3
  
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
    }
    else {
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
    
    # Clean this up:
    eig_out <- eig_replace( solve(S11) %*% S10 %*% solve(S00) %*% S01 )
    V <- eig_out$V
    D <- eig_out$D
    
    # V <- sortrows( [ V diag(D) ], p1+1 )
    V <- sortrows( whatever )
    V <- V[1:p1, ]
    # betaStar <- V( 1:p1, p1 : -1 : p1-r+1 )
    betaStar <- V[ 1:p1, seq(p1, p1-r+1, by = -1) ]
    
    
    
    
    # If either alpha or beta is restricted, then the likelihood is
    #   maximized subject to the constraints. This section of the code
    #   replaces the call to the switching algorithm in the previous
    #   version.
    if(!is.null(opt$R_Alpha) | !is.null(opt$R_Beta) ) {
      
      
      # [ betaStar, alphaHat, OmegaHat ]...
      # = RstrctOptm_Switch(betaStar, S00, S01, S11, T, p, opt)
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
      betaHat <- betaStar[1:p, 1:r]
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
  }
  else {
    if(r > 0) {
      GammaHat <- t( solve(t(Z2hat) %*% Z2hat) %*% t(Z2hat) %*% (Z0hat - Z1hat %*% t(PiStar)) )
    }
    else {
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



FCVARlikeMu <- function(y, db, mu, k, r, opt) {
  
  
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



FCVARlike <- function(x, params, k, r, opt) {
  
  # global estimatesTEMP
  
  T <- nrow(x)
  p <- ncol(x)
  
  
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
      mu <- params[2:end]
    }
    
  }
  else {
    d <- params[1]
    b <- params[2]    
    if (opt$levelParam) {
      mu <- params[3:end]
    }
    
  }
  
  
  # Modify the data by a mu level.
  if (opt$levelParam) {
    y <- x - matrix(1, nrow = T, ncol = p) %*% diag(mu)
  }
  else {
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
  }
  else {
    estimatesTEMP$muHat <- NULL
  }
  
  
  # Calculate value of likelihood function.
  T <- nrow(y) - opt$N
  p <- ncol(y)
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates$OmegaHat))
  
  
  return(like)
}







# Functions to make: 
# Z_array <- TransformData(x, k, db, opt)


#--------------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------------




################################################################################
# End
################################################################################

