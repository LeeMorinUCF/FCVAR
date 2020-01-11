################################################################################
# 
# Main Estimation Function for FCVAR
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
# DESCRIPTION: The main FCVAR estimation function. 
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
# Define function to estimate FCVAR model
################################################################################
# 
# function [ results ] <- FCVARestn(x,k,r,opt)
# Written by Michal Popiel and Morten Nielsen (This version 04.09.2016)
# 
# DESCRIPTION: This function performs estimation of the FCVAR system. It is
# 	the main function in the program with several nested functions, each
# 	described below. It estimates the model parameters, calculates the
# 	standard errors and the number of free parameters, obtains the residuals
# 	and the roots of the characteristic polynomial, and prints the output.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- results (a Matlab structure containing estimation results)
#            - results$startVals     (Starting values used for optimization)
#            - results$options       (Estimation options)
#            - results$like          (Model log-likelihood)
#            - results$coeffs        (Parameter estimates)
#            - results$rankJ         (Rank of Jacobian for
#					identification condition)
#            - results$fp            (Number of free parameters)
#            - results$SE            (Standard errors)
#            - results$NegInvHessian (Negative of inverse Hessian matrix)
#            - results$Residuals     (Model residuals)
#            - results$cPolyRoots    (Roots of characteristic polynomial)
# 
################################################################################

FCVARestn <- function(x,k,r,opt) {
  
  # global estimatesTEMP # Different in R.
  
  # --- Preliminary steps --- %
  T <- nrow(x) - opt$N # number of observations
  p <- ncol(x)         # number of variables
  
  # Update options based on initial user input.
  opt <- updateRestrictions(opt, p, r)
  
  # Clear previous instances of coefficient estimates. This is done
  #   because estimatesTemp is a global structure that will not be cleared
  #   automatically if estimation is interrupted.
  estimatesTEMP <- NULL
  
  
  #--------------------------------------------------------------------------------
  # GRID SEARCH
  #--------------------------------------------------------------------------------
  
  # Hide all warnings for grid search.
  # warning('off', 'all')
  
  # Perform grid search and store results as starting values for
  #   numerical optimization below.
  if(opt$gridSearch) {
    print(sprintf('\nRunning grid search over likelihood for k=%g, r=%g.\n',...
            k,r))
    print(sprintf('This computation can be slow.\n'))
    print(sprintf('Set opt$gridSearch <- 0 to skip it.\n'))
    opt$db0 <- LikeGrid(x, k, r, opt)
    # Change upper and lower bounds to limit the search to some small
    #   interval around the starting values.
    opt$UB_db[1:2] <- pmin(opt$db0[1:2] + c(0.1, 0.1), opt$dbMax)
    opt$LB_db[1:2] <- pmax(opt$db0[1:2] - c(0.1, 0.1), opt$dbMin)
  }
  else {
    # Call to GetBounds returns upper/lower bounds for (d,b) or
    #  depending on whether or not restrictions have been imposed.
    UB_LB_bounds <- GetBounds(obj)
    opt$UB_db <- UB_LB_bounds$UB
    opt$LB_db <- UB_LB_bounds$LB
    
  }
  
  
  
  # Turn warnings back on for main estimation.
  # warning('on','all')
  
  # Hide warnings related to algorithm changes in numerical optimization.
  # warning('off','optim:fminunc:SwitchingMethod')
  # warning('off','optimL:fminunc:SwitchingMethod')
  
  # Hide warnings related to singular matrix in numerical optimization.
  # warning('off','MATLAB:nearlySingularMatrix')
  
  
  
  #--------------------------------------------------------------------------------
  # ESTIMATION
  #--------------------------------------------------------------------------------
  
  
  # Store equality restrictions, inequality restrictions, upper and lower
  #   bounds, and starting values for estimation because they need to be
  #   adjusted based on the presence of level parameter or d=b restriction.
  Rpsi <- opt$R_psi
  rpsi <- opt$r_psi
  startVals <- opt$db0[1:2] # starting values for mu get added later.
  Cdb <- opt$C_db
  cdb <- opt$c_db
  
  # If Rpsi is empty, then optimization is over (d,b), otherwise it is
  #  over phi. If it is over phi, need to make adjustments to startVals
  #  and make Cdb and cdb empty.
  
  
  
  if(size(Rpsi,1)==1) {
    
    H_psi <- null(Rpsi)
    
    # Need to back out phi from given db0
    startVals <- solve(t(H_psi) %*% H_psi) %*% t(H_psi) %*% t(startVals[1:2])
    
    if(opt$gridSearch) {
      # Translate from d,b to phi.
      UB <- solve(t(H_psi) %*% H_psi) %*% t(H_psi) %*% t(opt$UB_db)
      LB <- solve(t(H_psi) %*% H_psi) %*% t(H_psi) %*% t(opt$LB_db)
    }
    else {
      # Otherwise GetBounds returns the values in terms of phi.
      UB <- opt$UB_db
      LB <- opt$LB_db
    }
    
    # Add warning about turning these off if non-empty?
    Cdb <- NULL
    cdb <- NULL
    
    # Turn off equality restrictions if only one restriction is imposed
    #  because they will be imposed inside the profile likelihood
    #  function that is being maximized.
    
    Rpsi <- NULL
    rpsi <- NULL
    
    
  } else {
    
    UB <- opt$UB_db
    LB <- opt$LB_db
    
  }
  
  
  
  
  # If estimation involves level parameter, make appropriate
  #   adjustments so that the dimensions of the restrictions correspond to
  #   the number of parameters being estimated, i.e. fill with zeros.
  if(opt$levelParam) {
    
    
    # If there are equality restrictions add coefficients of 0 to the
    #  level parameter.
    if(!is.null(Rpsi)) {
      # Check conformability:
      Rpsi <- rbind(Rpsi, matrix(0, nrow = nrow(Rpsi), ncol = p))
    }
    
    # If there are inequality restrictions add coefficients of 0 to the
    #  level parameter.
    if(!is.null(Cdb)) {
      # Check conformability:
      Cdb <- rbind(Cdb, matrix(0, nrow = nrow(Cdb), ncol = p))
    }
    
    # If grid search has not been used to set all initial values, add p
    #   starting values, set as the first observations on each variable,
    #   to the startVals variable to account for the level parameter.
    if(opt$gridSearch == 0) {
      # Check conformability:
      startVals <- c(startVals, x[1, ])
    }
    else {
      # Check conformability:
      startVals <- c(startVals, opt$db0[3:length(opt$db0)])
    }
    
    
    # Level parameter is unrestricted, but the length of UB and LB
    #   needs to match the number of parameters in optimization.
    UB <- c(UB, matrix(1,nrow = 1, ncol = p)*Inf)
    LB <- c(LB, -matrix(1,nrow = 1, ncol = p)*Inf)
    
    
  }
  
  
  
  if(size(opt$R_psi)==2) {
    
    
    # d,b are exactly identified by the linear restrictions and Rpsi is
    #  invertible. We use opt$R_psi here because Rpsi is adjusted
    #  depending on the presence of level parameters. Transpose is
    #  necessary to match input/output of other cases.
    dbTemp <- t(t(opt$R_psi) %*% solve(opt$R_psi*t(opt$R_psi)) %*% opt$r_psi)
    y <- x
    
    if(opt$levelParam) {
      # Optimize over level parameter for given (d,b).
      StartVal <- y[1, ]
      
      # [ muHat, maxLike, ! ] ...
      # <- fminunc(@( params ) -FCVARlikeMu(x, dbTemp, params, k, r, opt), ...
      #            StartVal, opt$UncFminOptions )
      
      # Need to implement optimization correctly. 
      min_out <- optim_unc(-FCVARlikeMu(params, x, dbTemp, k, r, opt), 
                           startVal, opt$UncFminOptions)
      
      muHat <- min_out$par
      maxLike <- min_out$f
      
      
      y <- x - matrix(1, nrow = T+opt$N, ncol = p) %*% diag(muHat)
    }
    else {
      maxLike <- -FCVARlike(dbTemp, y, k, r, opt)
    }
    
    
    # Obtain concentrated parameter estimates.
    estimates <- GetParams(y, k, r, dbTemp, opt)
    
    # Storing the estimates in a global structure here allows us to skip a
    #   call to GetParams after optimization to recover the coefficient
    #   estimates
    estimatesTEMP <- estimates
    # If level parameter is present, they are the last p parameters in the
    #   params vector
    if (opt$levelParam) {
      estimatesTEMP$muHat <- muHat
    }
    else {
      estimatesTEMP$muHat <- NULL
    }
    
    
  } else {
    # [ !, maxLike, ! ] ...
    # <- fmincon(@( params ) -FCVARlike(x, params, k, r, opt), ...
    #            startVals, Cdb, cdb, Rpsi, rpsi, LB, UB, [], opt$ConFminOptions )
    
    # Need to implement optimization correctly. 
    min_out <- optim_con(-FCVARlike(params, x, k, r, opt), 
                     startVals, Cdb, cdb, Rpsi, rpsi, LB, UB, opt$ConFminOptions)
    
    maxLike <- min_out$f
  }
  

  
  
  
  
  #--------------------------------------------------------------------------------
  # Store outputs
  #--------------------------------------------------------------------------------
  
  # Store the updated estimation options. 
  results$startVals <- startVals
  results$options <- opt
  results$options$UB_db <- UB
  results$options$LB_db <- LB
  
  
  # Adjust the sign of the likelihood and store the results
  maxLike <- -maxLike
  results$like <- maxLike
  
  # Coefficients are taken from a global defined in the likelihood
  #   function
  results$coeffs <- estimatesTEMP
  
  
  
  
  
  
  return(results)
}







      







################################################################################
# Define function to ...
################################################################################
# 

# 
################################################################################






################################################################################
# End
################################################################################
