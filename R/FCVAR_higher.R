################################################################################
# 
# Pre- and Post-Estimation Functions for FCVAR
# 
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 10, 2020
# 
################################################################################
# 
# Written by Michal Popiel and Morten Nielsen (This version 01.10.2020)
#
# DESCRIPTION: A library of functions that call the FCVAR estimation function. 
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
# Define function to select lag order.
################################################################################
# 
# function LagSelect(x, kmax, r, order, opt )
# Written by Michal Popiel and Morten Nielsen (This version 3.31.2016)
# 
# DESCRIPTION: This program takes a matrix of variables and performs lag
# 	selection on it by using the likelihood ratio test. Output and test
# 	results are printed to the screen.
#
# Input <- x     (matrix of variables to be included in the system)
#         kmax  (maximum number of lags)
#         r     (cointegration rank <- number of cointegrating vectors)
#         order (order of serial correlation for white noise tests)
#         opt   (object containing estimation options)
# Output <- none (only output to screen)
# 
################################################################################


LagSelect <- function(x, kmax, r, order, opt ) {
  
  
  # Determine (initial) dimensions of system.
  T <- nrow(x) # Length of sample (before truncation for initial values).
  p <- ncol(x) # Number of variables.
  printWN <- 0    # Set to zero to avoid printing output for each WN test
  
  # Do not print FCVAR estimation for each lag in the loop.
  opt$print2screen <- 0
  
  # Do not plot roots of characteristic polynomial for each lag in the loop.
  opt$plotRoots <- 0
  
  # Do not calculate standard errors.
  opt$CalcSE <- 0
  
  
  #--------------------------------------------------------------------------------
  # Create storage bins for output
  #--------------------------------------------------------------------------------
  
  # Row j of each of the following contains the associated results for
  # 	lag length j+1
  D        <- matrix(0, nrow = kmax+1, ncol = 2) # estimates of d and b
  loglik   <- matrix(0, nrow = kmax+1, ncol = 1) # log-likelihood
  LRtest   <- matrix(0, nrow = kmax+1, ncol = 1) # likelihood ratio test statistic for 
  #	significance of Gamma_{j+1}
  pvLRtest <- matrix(0, nrow = kmax+1, ncol = 1) # likelihood ratio test P-value
  aic      <- matrix(0, nrow = kmax+1, ncol = 1) # Akaike information criterion
  bic      <- matrix(0, nrow = kmax+1, ncol = 1) # Bayesian information criterion
  pvMVq    <- matrix(0, nrow = kmax+1, ncol = 1) # multivariate residual white noise Q-test P-value
  pvWNQ    <- matrix(0, nrow = kmax+1, ncol = p) # univariate residual white noise Q-test P-values
  #   for each residual
  pvWNLM   <- matrix(0, nrow = kmax+1, ncol = p) # univariate residual white noise LM-test P-values
  #   for each residual
  
  
  
  #--------------------------------------------------------------------------------
  # Estimate FCVAR for each lag order. 
  #--------------------------------------------------------------------------------
  
  
  for (k in 0:kmax) {
    
    # ----- Estimation ---------#         
    results <- FCVARestn(x,k,r,opt)
    
    # ----- Record relevant output ---------%
    loglik[k+1] <- results$like
    D[k+1, ]    <- results$coeffs.db
    aic[k+1]    <- -2*loglik[k+1] + 2*results$fp
    bic[k+1]    <- -2*loglik[k+1] + results$fp*log(T - opt$N)
    
    # ----- White noise tests ---------%
    pvWN_list <- mv_wntest(results$Residuals, order, printWN)   
    # [ ~, pvWNQ(k+1,:), ~, pvWNLM(k+1,:), ~, pvMVq(k+1,:) ]  <- mv_wntest(...)
    # Need to make sense of this after writing mv_wntest
    pvWNQ[k+1, ] <- pvWN_list[2]
    pvWNLM[k+1, ] <- pvWN_list[4]
    pvMVq[k+1, ] <- pvWN_list[6]
    
    # ----- LR test of lag <- k vs lag <- k-1 -----%
    if (k > 0) {
      LRtest(k+1)   <- 2*(loglik(k+1)-loglik(k))
      pvLRtest(k+1) <- 1-chi2cdf(LRtest(k+1), p^2)
    }
    
    
    
  }
  
  
  # Find lag corresponding to min of information criteria
  i_aic <- which.min(aic)
  i_bic <- which.min(bic)
  
  #--------------------------------------------------------------------------------
  # Print output 
  #--------------------------------------------------------------------------------
  
  # create a variable for output strings
  yesNo <- c('No','Yes') # Ironic order, No?
  
  print(sprintf('\n-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf('                        Lag Selection Results \n'))
  print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf('Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n', p, T))
  print(sprintf('Order for WN tests:   %6.0f     Number of observations for estimation:  %6.0f \n', order, T-opt$N))
  print(sprintf('Restricted constant:  %6s     Initial values:                         %6.0f\n', yesNo{opt$rConstant+1}, opt$N )   )
  print(sprintf('Unrestricted constant: %6s     Level parameter:                        %6s\n', yesNo{opt$unrConstant+1}, yesNo{opt$levelParam+1} ))
  print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf('k  r    d    b      LogL     LR    pv    AIC       BIC     pmvQ'))
  for (i in 1:p) {
    print(sprintf(' pQ%1.0f  pLM%1.0f', i,i))
  }
  
  
  print(sprintf('\n'))
  for (k in seq(kmax, 0, by = -1) ) {
    
    #       print(sprintf('%2.0f %2.0f %4.3f %4.3f %7.2f %6.2f %5.3f %8.2f %8.2f %4.2f',      
    #             k, r, D(k+1,:), loglik(k+1), LRtest(k+1),
    #           pvLRtest(k+1), aic(k+1), bic(k+1), pvMVq(k+1,:))
    print(sprintf('%2.0f %2.0f %4.3f %4.3f %7.2f %6.2f %5.3f %8.2f',
                  k, r, D[k+1, ], loglik[k+1], LRtest[k+1],
                  pvLRtest[k+1], aic[k+1]))
    # For AIC add asterisk if min value
    if(k+1 == i_aic) {print(sprintf('*'))} else {print(sprintf(' '))} 
    # Print BIC information criteria and add asterisk if min value
    print(sprintf(' %8.2f', bic[k+1]))
    if(k+1 == i_bic) {print(sprintf('*'))} else {print(sprintf(' '))}
    # Print multivariate white noise test P-values
    print(sprintf(' %4.2f', pvMVq[k+1, ]))
    # Print the individual series white noise test P-values
    for (i in 1:p) {
      print(sprintf(' %4.2f %4.2f', pvWNQ[k+1,i], pvWNLM[k+1,i]))
    }
    
    
    print(sprintf('\n'))
    
  }
  
  
  print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
  
}


################################################################################
# Define function to ...
################################################################################
# 
# function [ rankTestStats ] <- RankTests(x, k, opt)
# Written by Michal Popiel and Morten Nielsen (This version 11.17.2014)
# Based on Lee Morin & Morten Nielsen (June 5, 2013)
#
# DESCRIPTION: Performs a sequence of  likelihood ratio tests
# 	for cointegrating rank.
#
# The results are printed to screen if the indicator print2screen is 1.
#
# input <- vector or matrix x of data.
#       scalar k denoting lag length.
#       opt (object containing estimation options)
#
# output <- rankTestStats structure with results from cointegrating rank
#           tests, containing the following (p+1) vectors with ith element
#           corresponding to rank <- i-1:
#	dHat	(estimates of d)
#	bHat	(estimate of b)
#	LogL	(maximized log-likelihood)
#	LRstat  (LR trace statistic for testing rank r against rank p)
#	pv      (P-value of LR trace test, or "999" if P-value is not available)
# 
################################################################################

RankTests <- function(x, k, opt) {
  
  T <- nrow(x) - opt$N
  p <- ncol(x)
  
  # Store user specified options for printing to screen because it will be
  # turned off while looping over ranks.
  tempPrint2Screen <- opt$print2screen
  
  # Create output storage to be filled.
  bHat   <- matrix(0, nrow = p+1, ncol = 1)
  dHat   <- matrix(0, nrow = p+1, ncol = 1)
  LogL   <- matrix(0, nrow = p+1, ncol = 1)
  LRstat <- matrix(0, nrow = p+1, ncol = 1)
  pv     <- matrix(0, nrow = p+1, ncol = 1)
  
  # Do not print FCVAR estimation for each rank in the loop.
  opt$print2screen <- 0
  
  # Do not plot roots of characteristic polynomial for each lag in the loop.
  opt$plotRoots <- 0
  
  # Do not calculate standard errors.
  opt$CalcSE <- 0
  
  
  
  # For calculation of P-values
  if(opt$rConstant | opt$levelParam) {
    consT <- 1
  }
  else {
    consT <- 0
  }
  
  
  # Estimate models for all ranks
  for (r in 0 : p) {
    results <- FCVARestn(x,k,r, opt)
    dHat[r+1] <- results$coeffs.db[1]
    bHat[r+1] <- results$coeffs.db[2]
    LogL[r+1] <- results$like
  }
  
  
  
  
  # Calculate the LR statistics and P-values
  for (r in 0 : p-1) {
    
    LRstat[r+1] <-  - 2*( LogL[r+1] - LogL[p+1] )
    
    p_val = NULL
    # Get P-values, if
    # (1) no deterministic terms, or
    # (2) there is only restricted constant and d=b, or
    # (3) there is only a level parameter and d=b.
    if((~opt$rConstant & ~opt$unrConstant & ~opt$levelParam) |
       (opt$rConstant  & ~opt$unrConstant & opt$restrictDB) |
       (opt$levelParam & ~opt$unrConstant & opt$restrictDB) ) {
      p_val <- get_pvalues(p-r, bHat[r+1], consT, LRstat[r+1], opt)
    }
    
    # If automatic calls to P-value program have not been installed or
    # enabled, then p_val is empty. Set it to 999 so that it can have a
    # value for storage in the rankTestStats matrix below.
    if(is.null(p_val)) {
      p_val <- 999
    }
    
    # Store P-values.
    pv[r+1] <- p_val
    
  }
  
  
  
  
  # Restore settings.
  opt$print2screen <- tempPrint2Screen
  
  # Print the results to screen.
  if opt$print2screen {
    
    # create a variable for output strings
    yesNo <- c('No','Yes')
    
    print(sprintf('\n-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('                         Likelihood Ratio Tests for Cointegrating Rank                               \n'))
    print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n', p, T+opt$N))
    print(sprintf('Number of lags:       %6.0f     Number of observations for estimation:  %6.0f \n', k, T))
    print(sprintf('Restricted constant:  %6s     Initial values:                         %6.0f\n', yesNo{opt$rConstant+1}, opt$N ))
    print(sprintf('Unestricted constant: %6s     Level parameter:                        %6s\n', yesNo{opt$unrConstant+1}, yesNo{opt$levelParam+1} ))
    print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('Rank \t  d  \t  b  \t Log-likelihood\t LR statistic\t P-value\n'))
    for (i in 1:p) {
      if (pv(i) != 999) {
        print(sprintf('%2.0f   \t%5.3f\t%5.3f\t%15.3f\t%13.3f\t%8.3f\n', i-1, dHat(i), bHat(i), LogL(i), LRstat(i), pv(i)))
      }
      else {
        print(sprintf('%2.0f   \t%5.3f\t%5.3f\t%15.3f\t%13.3f\t    ----\n', i-1, dHat(i), bHat(i), LogL(i), LRstat(i)))
      }
      
    }
    
    print(sprintf('%2.0f   \t%5.3f\t%5.3f\t%15.3f\t         ----\t    ----\n', i, dHat(i+1), bHat(i+1), LogL(i+1)))
    print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
    
    
  }
  
  
  
  
  
  
  # Return structure of rank test results$
  rankTestStats$dHat   <- dHat
  rankTestStats$bHat   <- bHat
  rankTestStats$LogL   <- LogL
  rankTestStats$LRstat <- LRstat
  rankTestStats$pv     <- pv
  
  
  return(rankTestStats)
}
end



################################################################################
# Define function to ...
################################################################################
# 
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
#
# DESCRIPTION: This function calls the program FDPVAL in the terminal and
# 	returns the P-value based on the user's inputs. The function's
# 	arguments must be converted to strings in order to interact with the
# 	terminal.
#
# Input <- q        (number of variables minus rank)
#         b        (parameter)
#         consT    (boolean variable indicating whether or not there is
#                   constant present)
#         testStat (value of the test statistic)
#	      opt (object containing estimation options)
# Output <- pv (P-value for likelihood ratio test)
# 
################################################################################


get_pvalues <- function(q, b, consT, testStat, opt) {
  
  
  
  if(b < 0.5) {
    # Series are stationary so use chi^2 with (p-r)^2 df, see JN2012
    pv <- 1-chi2cdf(testStat, q^2)
  } else {
    
    # For non-stationary systems use simulated P-values from
    # MacKinnon and Nielsen (2014, Journal of Applied Econometrics)
    # and the C++ program conversion by Jason Rhinelander.
    
    # Convert user input to system commands and
    # translate arguments to strings
    args <- sprintf('%g %g %g %g', q, b, consT, testStat)
    # Combine path, program, and arguments
    outCode <- sprintf('%s %s', opt$progLoc, args)
    
    # Evaluate system command
    # Note: fdpval is a separately installed program.
    # For more information see: https://github.com/jagerman/fracdist
    # For download see https://github.com/jagerman/fracdist/releases
    [ flag , pval] <- system([outCode])
    # Doesn't work on all platforms. 
    # Should replace with a wrapper for the fracdist code. 
    
    # Note: string is returned, so it needs to be converted
    if(flag==0) {
      
      # check if program was executed without errors
      pv <- str2double(pval)
      
    } else {
      
      # program failed
      pv <- NULL
      
    } 
    
  }
  
  return(pv)
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
