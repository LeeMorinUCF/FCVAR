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

  # Do not print output for each WN test.
  printWN <- 0

  # Do not print FCVAR estimation for each lag in the loop.
  opt$print2screen <- 0

  # Do not plot roots of characteristic polynomial for each lag in the loop.
  opt$plotRoots <- 0

  # Do not calculate standard errors.
  opt$CalcSE <- 0

  # Shouldn't this be done here, too?
  # Update options based on initial user input.
  # opt <- updateRestrictions(opt, p, r)
  # Maybe not, because it is already done in FCVARestn.


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

    cat(sprintf('Estimating for k = %d and r = %d.\n\n', k, r))

    # ----- Estimation ---------#
    results <- FCVARestn(x, k, r, opt)


    cat(sprintf('Finished Estimation for k = %d and r = %d.\n\n', k, r))

    # ----- Record relevant output ---------%
    loglik[k+1] <- results$like
    D[k+1, ]    <- results$coeffs$db
    aic[k+1]    <- -2*loglik[k+1] + 2*results$fp
    bic[k+1]    <- -2*loglik[k+1] + results$fp*log(T - opt$N)

    # ----- White noise tests ---------%
    mv_wntest_out <- mv_wntest(results$Residuals, order, printWN)
    # [ ~, pvWNQ(k+1,:), ~, pvWNLM(k+1,:), ~, pvMVq(k+1,:) ]  <- mv_wntest(...)
    # Need to make sense of this after writing mv_wntest. Check.
    pvWNQ[k+1, ] <- mv_wntest_out$pvQ
    pvWNLM[k+1, ] <- mv_wntest_out$pvLM
    pvMVq[k+1, ] <- mv_wntest_out$pvMVQ


    # ----- LR test of lag <- k vs lag <- k-1 -----%
    if (k > 0) {
      LRtest[k+1]   <- 2*(loglik[k+1] - loglik[k])
      pvLRtest[k+1] <- 1 - pchisq(LRtest[k+1], p^2)
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

  cat(sprintf('\n-----------------------------------------------------------------------------------------------------\n'))
  cat(sprintf('                        Lag Selection Results \n'))
  cat(sprintf('-----------------------------------------------------------------------------------------------------\n'))
  cat(sprintf('Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n', p, T))
  cat(sprintf('Order for WN tests:   %6.0f     Number of observations for estimation:  %6.0f \n', order, T-opt$N))
  cat(sprintf('Restricted constant:  %6s     Initial values:                         %6.0f\n', yesNo[opt$rConstant+1], opt$N )   )
  cat(sprintf('Unrestricted constant: %6s     Level parameter:                        %6s\n', yesNo[opt$unrConstant+1], yesNo[opt$levelParam+1] ))
  cat(sprintf('-----------------------------------------------------------------------------------------------------\n'))
  cat(sprintf('k  r    d    b      LogL     LR    pv    AIC       BIC     pmvQ'))
  for (i in 1:p) {
    cat(sprintf(' pQ%1.0f  pLM%1.0f', i,i))
  }


  cat(sprintf('\n'))
  for (k in seq(kmax, 0, by = -1) ) {

    #       cat(sprintf('%2.0f %2.0f %4.3f %4.3f %7.2f %6.2f %5.3f %8.2f %8.2f %4.2f',
    #             k, r, D(k+1,:), loglik(k+1), LRtest(k+1),
    #           pvLRtest(k+1), aic(k+1), bic(k+1), pvMVq(k+1,:))
    cat(sprintf('%2.0f %2.0f %4.3f %4.3f %7.2f %6.2f %5.3f %8.2f',
                  k, r, D[k+1, 1], D[k+1, 2], loglik[k+1], LRtest[k+1],
                  pvLRtest[k+1], aic[k+1]))
    # For AIC add asterisk if min value
    if(k+1 == i_aic) {cat(sprintf('*'))} else {cat(sprintf(' '))}
    # Print BIC information criteria and add asterisk if min value
    cat(sprintf(' %8.2f', bic[k+1]))
    if(k+1 == i_bic) {cat(sprintf('*'))} else {cat(sprintf(' '))}
    # Print multivariate white noise test P-values
    cat(sprintf(' %4.2f', pvMVq[k+1, ]))
    # Print the individual series white noise test P-values
    for (i in 1:p) {
      cat(sprintf(' %4.2f %4.2f', pvWNQ[k+1,i], pvWNLM[k+1,i]))
    }


    cat(sprintf('\n'))

  }


  cat(sprintf('-----------------------------------------------------------------------------------------------------\n'))

}


################################################################################
# Define function to perform a sequence of likelihood ratio tests
# 	for cointegrating rank.
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
  } else {
    consT <- 0
  }


  # Estimate models for all ranks
  for (r in 0 : p) {


    cat(sprintf('Estimating for k = %d and r = %d.\n\n', k, r))

    results <- FCVARestn(x, k, r, opt)

    cat(sprintf('Finished Estimation for k = %d and r = %d.\n\n', k, r))

    dHat[r+1] <- results$coeffs$db[1]
    bHat[r+1] <- results$coeffs$db[2]
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
    if (FALSE & (
      (!opt$rConstant & !opt$unrConstant & !opt$levelParam) |
      (opt$rConstant  & !opt$unrConstant & opt$restrictDB) |
      (opt$levelParam & !opt$unrConstant & opt$restrictDB) )  ) {

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
  if (opt$print2screen) {

    # create a variable for output strings
    yesNo <- c('No','Yes')

    cat(sprintf('\n-----------------------------------------------------------------------------------------------------\n'))
    cat(sprintf('                         Likelihood Ratio Tests for Cointegrating Rank                               \n'))
    cat(sprintf('-----------------------------------------------------------------------------------------------------\n'))
    cat(sprintf('Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n', p, T+opt$N))
    cat(sprintf('Number of lags:       %6.0f     Number of observations for estimation:  %6.0f \n', k, T))
    cat(sprintf('Restricted constant:  %6s     Initial values:                         %6.0f\n', yesNo[opt$rConstant+1], opt$N ))
    cat(sprintf('Unestricted constant: %6s     Level parameter:                        %6s\n', yesNo[opt$unrConstant+1], yesNo[opt$levelParam+1] ))
    cat(sprintf('-----------------------------------------------------------------------------------------------------\n'))
    cat(sprintf('Rank \t  d  \t  b  \t Log-likelihood\t LR statistic\t P-value\n'))
    for (i in 1:p) {
      if (pv[i] != 999) {
        cat(sprintf('%2.0f   \t%5.3f\t%5.3f\t%15.3f\t%13.3f\t%8.3f\n', i-1, dHat[i], bHat[i], LogL[i], LRstat[i], pv[i]))
      }
      else {
        cat(sprintf('%2.0f   \t%5.3f\t%5.3f\t%15.3f\t%13.3f\t    ----\n', i-1, dHat[i], bHat[i], LogL[i], LRstat[i]))
      }

    }

    cat(sprintf('%2.0f   \t%5.3f\t%5.3f\t%15.3f\t         ----\t    ----\n', i, dHat[i+1], bHat[i+1], LogL[i+1]))
    cat(sprintf('-----------------------------------------------------------------------------------------------------\n'))


  }






  # Return list of rank test results.
  rankTestStats <- list(
    dHat   = dHat,
    bHat   = bHat,
    LogL   = LogL,
    LRstat = LRstat,
    pv     = pv
  )



  return(rankTestStats)
}




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

  print('b = ')
  print(b)

  if (b < 0.5) {
    # Series are stationary so use chi^2 with (p-r)^2 df, see JN2012
    pv <- 1 - chi2cdf(testStat, q^2)
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
    # [ flag , pval] <- system([outCode])
    fracdist_out <- 7
    # Doesn't work on all platforms.
    # Should replace with a wrapper for the fracdist code.

    # Note: string is returned, so it needs to be converted
    if(flag == 0) {

      # check if program was executed without errors
      # pv <- str2double(pval)
      pv <- pval

    } else {

      # program failed
      pv <- NULL

    }

  }

  return(pv)
}
















################################################################################
# Define function to perform a multivariate Ljung-Box Q-test for
# 	white noise.
################################################################################
#
# function [ Q, pvQ, LM, pvLM, mvQ, pvMVQ ] <-
#                       mv_wntest(x, maxlag, printResults)
# Written by Michal Popiel and Morten Nielsen (This version 7.21.2015)
#
# DESCRIPTION: This function performs a multivariate Ljung-Box Q-test for
# 	white noise and univariate Q-tests and LM-tests for white noise on the
# 	columns of x.
# 	The LM test should be consistent for heteroskedastic series, Q-test is not.
#
# Input <- x            (matrix of variables under test, typically model residuals)
#         maxlag       (number of lags for serial correlation tests)
#         printResults (set =1 to print results to screen)
# Output <- Q     (1xp vector of Q statistics for individual series)
#		       pvQ 	 (1xp vector of P-values for Q-test on individual series)
#          LM    (1xp vector of LM statistics for individual series)
#		       pvLM	 (1xp vector of P-values for LM-test on individual series)
#          mvQ   (multivariate Q statistic)
#          pvMVQ (P-value for multivariate Q-statistic using p^2*maxlag df)
#
################################################################################


mv_wntest <- function(x, maxlag, printResults) {

  T <- nrow(x)
  p <- ncol(x)

  # Create bins for values
  pvQ  <- matrix(1, nrow = 1, ncol = p)
  pvLM <- matrix(1, nrow = 1, ncol = p)
  Q  <- matrix(0, nrow = 1, ncol = p)
  LM <- matrix(0, nrow = 1, ncol = p)

  # Perform univariate Q and LM tests and store the results.
  for (i in 1:p) {

    # Qtest_out <- Qtest(matrix(x[,i], nrow = T, ncol = 1), maxlag)
    Qtest_out <- Qtest(x[,i, drop = FALSE], maxlag)
    Q[i] <- Qtest_out$Qstat
    pvQ[i] <- Qtest_out$pv

    # LMtest_out <- LMtest(matrix(x[,i], nrow = T, ncol = 1), maxlag)
    LMtest_out <- LMtest(x[,i, drop = FALSE], maxlag)
    LM[i] <- LMtest_out$LMstat
    pvLM[i] <- LMtest_out$pv

  }

  # print('Finished the loop on columns')

  # Perform multivariate Q test.
  # [mvQ, pvMVQ] <- Qtest(x,maxlag)
  # Qtest_out <- Qtest(matrix(x, nrow = T, ncol = p), maxlag)
  Qtest_out <- Qtest(x[ , , drop = FALSE], maxlag)
  mvQ <- Qtest_out$Qstat
  pvMVQ <- Qtest_out$pv


  # Print output
  if (printResults) {
    cat(sprintf('\n       White Noise Test Results (lag = %g)\n', maxlag))
    cat(sprintf('---------------------------------------------\n'))
    cat(sprintf('Variable |       Q  P-val |      LM  P-val  |\n'))
    cat(sprintf('---------------------------------------------\n'))
    cat(sprintf('Multivar | %7.3f  %4.3f |     ----  ----  |\n', mvQ, pvMVQ))
    for (i in 1:p) {
      cat(sprintf('Var%g     | %7.3f  %4.3f | %7.3f  %4.3f  |\n',
              i, Q[i], pvQ[i], LM[i], pvLM[i] ))
    }

    cat(sprintf('---------------------------------------------\n'))
  }


  # Output a list of results.
  mv_wntest_out <- list(
    Q = Q,
    pvQ = pvQ,
    LM = LM,
    pvLM = pvLM,
    mvQ = mvQ,
    pvMVQ = pvMVQ
  )

  return(mv_wntest_out)
}



################################################################################
# Define function to perform a Breusch-Godfrey Lagrange Multiplier
# test for serial correlation
################################################################################
#
# Breusch-Godfrey Lagrange Multiplier test for serial correlation.
#
################################################################################


LMtest <- function(x,q) {

  # print('In LMtest')
  # print(size(x))
  # print(q)
  # print(nrow(x))
  # print(ncol(x))
  # print('(q+1):T = ')
  # print((q+1):T)


  # cat(sprintf('x is %d by %d', nrow(x), ncol(x)))

  # Breusch-Godfrey Lagrange Multiplier test for serial correlation.
  T <- nrow(x)
  x <- x - mean(x)

  # cat(sprintf('x is %d by %d', nrow(x), ncol(x)))

  # print(x)

  y <- x[seq(q+1, T), ,  drop = FALSE]
  z <- x[seq(1, T-q), ,  drop = FALSE]

  for (i in 1:(q-1)) {
    z <- cbind(x[seq(i+1, T-q+i), ,  drop = FALSE], z)
  }



  # print(head(y))
  # print(tail(y))
  # cat(sprintf('y is %d by %d', nrow(y), ncol(y)))


  e <- y
  # print(head(e))
  # print(tail(e))


  # cat(sprintf('e is %d by %d', nrow(e), ncol(e)))
  # cat(sprintf('z is %d by %d', nrow(z), ncol(z)))
  # kron_test <- kron(matrix(1, 1, q), e)
  # cat(sprintf('kron_test is %d by %d', nrow(kron_test), ncol(kron_test)))


  # s <- z[,1:q] * repmat(e,1,q)
  # Translate this properly:
  s <- z[,1:q,  drop = FALSE] * kron(matrix(1, 1, q), e)

  # cat(sprintf('s is %d by %d', nrow(s), ncol(s)))

  sbar <- t(colMeans(s))

  # print('sbar = ')
  # print(sbar)

  kron_sbar <- kron(matrix(1, nrow(s)), sbar)

  # print(head(kron_sbar))
  # print(tail(kron_sbar))

  # The next line bsxfun(@FUNC, A, B) applies the element-by-element binary
  # operation FUNC to arrays A and B, with singleton expansion enabled.
  # Need to translate this to R:
  # s <- mapply(minus, s, sbar)
  # Never mind that fancy pants stuff for something that
  # is not even calculated in a loop.
  s <- s - kron_sbar

  S <- t(s) %*% s/T
  # LMstat <- T*sbar %*% S^(-1) %*% t(sbar)
  LMstat <- T*sbar %*% solve(S) %*% t(sbar)
  pv <- 1 - pchisq(LMstat, q)


  # Output a list of results.
  LMtest_out <- list(
    LMstat = LMstat,
    pv = pv
  )

  return(LMtest_out)
}



################################################################################
# Define function to perform a Ljung-Box Q-test for serial correlation
################################################################################
#
# (Multivariate) Ljung-Box Q-test for serial correlation, see
# 	Luetkepohl (2005, New Introduction to Multiple Time Series Analysis, p. 169).
#
################################################################################

Qtest <- function(x, maxlag) {

  # print('class(x) = ')
  # print(class(x))

  T <- nrow(x)
  p <- ncol(x)

  # print('p = ')
  # print(p)
  # print('T = ')
  # print(T)

  # t <- 7
  # print(x[t, , drop = FALSE])
  # print(t(x[t, , drop = FALSE]))
  # print(t(x[t, , drop = FALSE]) %*% x[t, , drop = FALSE])

  C0 <- matrix(0, nrow = p, ncol = p)
  for (t in 1:T) {
    # C0 <- C0 + x[t, ] %*% t(x[t, ])
    C0 <- C0 + t(x[t, , drop = FALSE]) %*% x[t, , drop = FALSE]
  }
  C0 <- C0/T

  C <- array(rep(0, p*p*maxlag), dim = c(p,p,maxlag))
  for (i in 1:maxlag) {
    for (t in (i+1):T) {
      C[ , ,i] <- C[ , ,i] + t(x[t, , drop = FALSE]) %*% x[t-i, , drop = FALSE]
    }

    C[ , ,i] <- C[ , ,i]/(T-i) # Note division by (T-i) instead of T.
  }


  # (Multivariate) Q statistic
  Qstat <- 0
  for (j in 1:maxlag) {
    # Qstat <- Qstat+trace(C(:,:,j)'*inv(C0)*C(:,:,j)*inv(C0)) / (T-j) %'
    # The following line is a more efficient calculation than the previous
    # Qstat <- Qstat+trace( (C(:,:,j)'/C0)*(C(:,:,j)/C0) ) / (T-j) %' #'
    # Need function for trace:
    Qstat <- Qstat + sum(diag( (t(C[ , ,j]) %*% solve(C0)) %*% (C[ , ,j] %*% solve(C0)) )) / (T-j)
  }


  Qstat <- Qstat*T*(T+2)
  pv <- 1 - pchisq(Qstat, p*p*maxlag) # P-value is calculated with p^2*maxlag df.


  # Output a list of results.
  Qtest_out <- list(
    Qstat = Qstat,
    pv = pv
  )

  return(Qtest_out)
}





################################################################################
# Define function to perform a likelihood ratio test
################################################################################
#
# function results <- HypoTest(modelUNR, modelR)
# Written by Michal Popiel and Morten Nielsen (This version 2.24.2015)
#
# DESCRIPTION: This function performs a likelihood ratio test of the null
# 	hypothesis: "model is modelR" against the alternative hypothesis:
# 	"model is modelUNR".
#
# Input <- modelUNR (structure of estimation results created for unrestricted model)
#         modelR (structure of estimation results created for restricted model)
# Output <- results: a Matlab structure containing test results
#            - results$loglikUNR (loglikelihood of unrestricted model)
#            - results$loglikR   (loglikelihood of restricted model)
#            - results$df        (degrees of freedom for the test)
#            - results$LRstat    (likelihood ratio test statistic)
#            - results$p_LRtest  (P-value for test)
#
################################################################################


HypoTest <- function(modelUNR, modelR) {


  # Calculate the test statistic.
  LR_test <- 2*(modelUNR$like - modelR$like)

  # Calculate the degrees of freedom by taking the difference in free
  # parameters between the unrestricted and restricted model.
  df <- modelUNR$fp - modelR$fp

  # Calculate the P-value for the test.
  p_LRtest <- 1 - pchisq(LR_test, df)

  # Print output.
  cat(sprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n',
                modelUNR$like, modelR$like))
  cat(sprintf('Test results (df <- %1.0f):\nLR statistic: \t %3.3f\nP-value: \t %1.3f\n',
                df,LR_test,p_LRtest))


  # Return the test results in a list.
  LRtest <- list(
    loglikUNR = modelUNR$like,
    loglikR   = modelR$like,
    df        = df,
    LRstat    = LR_test,
    pv        = p_LRtest
  )

  return(LRtest)
}


################################################################################
# Define function to simulate the FCVAR model
################################################################################
#
# function xSim <- FCVARsim(data, model, NumPeriods)
# Written by Michal Popiel and Morten Nielsen (This version 08.06.2015)
#
# DESCRIPTION: This function simulates the FCVAR model as specified by
#               input "model" and starting values specified by "data."
#               Errors are drawn from a Normal distribution.
#
# Input <- data       (N x p matrix of data)
#         model      (a Matlab structure containing estimation results)
#         NumPeriods (number of steps for simulation)
# Output <- xSim      (NumPeriods x p matrix of simulated data values)
#
################################################################################

FCVARsim <- function(x, model, NumPeriods) {


  #--------------------------------------------------------------------------------
  # Preliminary definitions
  #--------------------------------------------------------------------------------

  # x <- data
  # p <- ncol(data)
  p <- ncol(x)
  opt <- model$options
  cf  <- model$coeffs
  d <- cf$db[1]
  b <- cf$db[2]

  # Generate disturbance term
  err <- matrix(rnorm(NumPeriods*p), nrow = NumPeriods, ncol = p)

  #--------------------------------------------------------------------------------
  # Recursively generate simulated data values
  #--------------------------------------------------------------------------------

  for (i in 1:NumPeriods) {


    # append x with zeros to simplify calculations.
    x <- rbind(x, rep(0, p))
    T <- nrow(x)

    # Adjust by level parameter if present.
    if(opt$levelParam) {
      y <- x - matrix(1, nrow = T, ncol = 1) %*% cf$muHat
    } else {
      y <- x
    }


    # Main term, take fractional lag.
    z <- Lbk(y,d,1)

    # Error correction term.
    if (!is.null(cf$alphaHat)) {

      z <- z + FracDiff( Lbk(y, b, 1), d - b ) %*% t(cf$PiHat)

            if (opt$rConstant) {
              z <- z + FracDiff( Lbk(matrix(1, nrow = T, ncol = 1), b, 1), d - b ) %*%
                cf$rhoHat %*% t(cf$alphaHat)
            }

    }


    # Add unrestricted constant if present.
    if(opt$unrConstant) {
      z <- z + matrix(1, nrow = T, ncol = 1) %*% t(cf$xiHat)
    }


    # Add lags if present.
    if (!is.null(cf$GammaHat)) {
      k <- ncol(cf$GammaHat)/p
      z <- z +  FracDiff(  Lbk( y , b, k)  , d) %*% t(cf$GammaHat)
    }


    # Adjust by level parameter if present.
    if (opt$levelParam) {
      z <- z + matrix(1, nrow = T, ncol = 1) %*% cf$muHat
    }


    # Add disturbance term
    z[T,] <- z[T,] + err[i, ]

    # Append to x matrix.
    x <- rbind(x[1:(T-1), ], z[T,])

  }

  #--------------------------------------------------------------------------------
  # Return simulated data values, including (EXcluding?) initial values.
  #--------------------------------------------------------------------------------

  xSim <- x[(nrow(data)+1):nrow(x), ]

  return(xSim)
}


################################################################################
# Define function to simulate bootstrap samples for the FCVAR model
################################################################################
#
# function xBS <- FCVARsimBS(data, model, NumPeriods)
# Written by Michal Popiel and Morten Nielsen (This version 08.06.2015)
#
# DESCRIPTION: This function simulates the FCVAR model as specified by
#               input "model" and starting values specified by "data." It
#               creates a bootstrap sample by augmenting each iteration
#               with a bootstrap error. The errors are sampled from the
#               residuals specified under the "model" input and have a
#               positive or negative sign with equal probability
#               (Rademacher distribution).
#
# Input <- data       (T x p matrix of data)
#         model      (a Matlab structure containing estimation results)
#         NumPeriods (number of steps for simulation)
# Output <- xBS       (NumPeriods x p matrix of simulated bootstrap values)
#
################################################################################



FCVARsimBS <- function(data, model, NumPeriods) {


  #--------------------------------------------------------------------------------
  # Preliminary definitions
  #--------------------------------------------------------------------------------

  x <- data
  p <- ncol(data)
  opt <- model$options
  cf  <- model$coeffs
  d <- cf$db[1]
  b <- cf$db[2]
  T <- nrow(model$Residuals)

  # print('model = ')
  # print(model)

  # Generate disturbance term for use in the bootstrap.

  # print('ncol(model$Residuals) = ')
  # print(ncol(model$Residuals))
  # print('mean(model$Residuals) = ')
  # print(mean(model$Residuals))
  # print('colMeans(model$Residuals) = ')
  # print(colMeans(model$Residuals))
  # Duh!

  # Centre residuals
  res <- model$Residuals -
    matrix(1, nrow = nrow(model$Residuals), ncol = 1) %*% colMeans(model$Residuals)

  # Generate draws from Rademacher distribution for Wild bootstrap
  eRD <- - matrix(1, nrow = T, ncol = p) +
    2*( matrix(rnorm(T), nrow = T, ncol = 1) > 0) %*% matrix(1, nrow = 1, ncol = p)

  # Generate error term
  err <- res * eRD



  #--------------------------------------------------------------------------------
  # Recursively generate bootstrap sample
  #--------------------------------------------------------------------------------


  for (i in 1:NumPeriods) {


    # append x with zeros to simplify calculations
    x <- rbind(x, rep(0, p))
    T <- nrow(x)

    # Adjust by level parameter if present
    if(opt$levelParam) {
      y <- x - matrix(1, nrow = T, ncol = 1) %*% cf$muHat
    } else {
      y <- x
    }


    # Main term, take fractional lag
    z <- Lbk(y, d, 1)

    # print('cf$PiHat = ')
    # print(cf$PiHat)

    # print('!is.null(cf$alphaHat) = ')
    # print(!is.null(cf$alphaHat))
    # print('!is.na(cf$alphaHat) = ')
    # print(!is.na(cf$alphaHat))
    # print('!is.na(cf$alphaHat[1]) = ')
    # print(!is.na(cf$alphaHat[1]))

    # Error correction term
    if( !is.null(cf$alphaHat) && !is.na(cf$alphaHat[1])) {
      z <- z + FracDiff( Lbk(y, b, 1), d - b ) %*% t(cf$PiHat)
      if(opt$rConstant) {
        z <- z + FracDiff( Lbk(matrix(1, nrow = T, ncol = 1), b, 1), d - b ) %*%
          cf$rhoHat %*% t(cf$alphaHat)
      }

    }


    # Add unrestricted constant if present
    if(opt$unrConstant) {
      z <- z + matrix(1, nrow = T, nol = 1) %*% t(cf$xiHat)
    }


    # Add lags if present
    if(!is.null(cf$GammaHat)) {
      k <- ncol(cf$GammaHat)/p
      z <- z +  FracDiff(  Lbk( y , b, k)  , d) %*% t(cf$GammaHat)
    }


    # Adjust by level parameter if present
    if(opt$levelParam) {
      z <- z + matrix(1, nrow = T,ncol = 1) %*% cf$muHat
    }


    # Add disturbance term
    z[T, ] <- z[T, ] + err[i, ]

    # Append generated observation to x matrix
    x <- rbind(x[1:(T-1), ], z[T, ])

  }


  #--------------------------------------------------------------------------------
  # Return bootstrap sample
  #--------------------------------------------------------------------------------

  xBS <- x[(nrow(data)+1):nrow(x), ]

  return(xBS)
}


################################################################################
# Define function to calculate recursive forecasts with the FCVAR model
################################################################################
#
# function xf <- FCVARforecast(data, model, NumPeriods)
# Written by Michal Popiel and Morten Nielsen (This version 11.17.2014)
#
# DESCRIPTION: This function calculates recursive forecasts. It uses
#   FracDiff() and Lbk(), which are nested below.
#
# Input <- data (T x p matrix of data)
#         model (a Matlab structure containing estimation results)
#         NumPeriods (number of steps ahead for forecast)
# Output <- xf (NumPeriods x p matrix of forecasted values)
#
################################################################################



FCVARforecast <- function(x, model, NumPeriods) {


  #--------------------------------------------------------------------------------
  # Preliminary steps
  #--------------------------------------------------------------------------------

  # x <- x1
  # x <- data
  # p <- ncol(data)
  p <- ncol(x)
  opt <- model$options
  cf  <- model$coeffs
  d <- cf$db[1]
  b <- cf$db[2]

  #--------------------------------------------------------------------------------
  # Recursively generate forecasts
  #--------------------------------------------------------------------------------

  for (i in 1:NumPeriods) {


    # Append x with zeros to simplify calculations.
    # x <- rbind(x, matrix(0, nrow = 1, ncol = p))
    x <- rbind(x, rep(0, p))
    T <- nrow(x)

    # Adjust by level parameter if present.
    if(opt$levelParam) {
      y <- x - matrix(1, nrow = T, ncol = 1) %*% cf$muHat
    } else {
      y <- x
    }


    # Main term, take fractional lag.
    z <- Lbk(y, d, 1)

    # Error correction term.
    if(!is.null(cf$alphaHat)) {

      # print('size(cf$PiHat) = ')
      # print(size(cf$PiHat))
      # print('size(FracDiff( Lbk(y, b, 1), d - b )) = ')
      # print(size(FracDiff( Lbk(y, b, 1), d - b )))


      z <- z + FracDiff( Lbk(y, b, 1), d - b ) %*% t(cf$PiHat)
      if(opt$rConstant) {
        z <- z + FracDiff( Lbk(matrix(1, nrow = T, ncol = 1), b, 1), d - b ) %*%
          cf$rhoHat %*% t(cf$alphaHat)
      }

    }


    # Add unrestricted constant if present.
    if(opt$unrConstant) {
      z <- z + matrix(1, nrow = T, ncol = 1) %*% t(cf$xiHat)
    }


    # Add lags if present.
    if(!is.null(cf$GammaHat)) {
      k <- size(cf$GammaHat,2)/p
      z <- z +  FracDiff(  Lbk( y , b, k)  , d) %*% t(cf$GammaHat)
    }


    # Adjust by level parameter if present.
    if(opt$levelParam) {
      z <- z + matrix(1, nrow = T, ncol = 1) %*% cf$muHat
    }


    # Append forecast to x matrix.
    x <- rbind(x[1:(T-1),], z[T, ])

  }


  #--------------------------------------------------------------------------------
  # Return forecasts.
  #--------------------------------------------------------------------------------

  xf <- x[(nrow(data)+1):nrow(x), ]

  return(xf)
}


################################################################################
# Define function to calculate bootstrap likelihood ratio test statistics
################################################################################
#
# function [LRbs, H, mBS, mUNR] <- FCVARboot(x, k, r, optRES, optUNR, B)
# Written by Michal Popiel and Morten Nielsen (This version 08.06.2015)
#
# DESCRIPTION: This function generates a distribution of a likelihood ratio
#           test statistic using a Wild bootstrap, following the method of
#			      Boswijk, Cavaliere, Rahbek, and Taylor (2013). It takes two sets
#           of options as inputs to estimate the model under the null and the
#           unrestricted model.
#
# Input <- x      (data - if k>0, actual data is used for initial values)
#         k      (number of lags)
#         optRES (options object for restricted model under the null)
#         optUNR (options object to estimate unrestricted model)
#         B      (number of bootstrap samples)
#
# Output <- LRbs (B x 1 vector simulated likelihood ratio statistics)
#          pv   (approximate p-value for LRstat based on bootstrap
#                                                       distribution)
#          H    (a Matlab structure containing LR test results, it is
#               identical to the output from HypoTest, with one addition,
#               namely H$pvBS which is the Bootstrap P-value)
#          mBS  (model estimates under the null)
#          mUNR (model estimates under the alternative)
#
################################################################################


FCVARboot <- function(x, k, r, optRES, optUNR, B) {


  # Calculate length of sample to generate, adjusting for initial values
  T <- nrow(x) - optRES$N

  # Use first k+1 observations for initial values
  data <- x[1:(k+1), ]

  LR <- matrix(0, nrow = B, ncol = 1)

  # Turn off output and calculation of standard errors for faster computation
  optUNR$print2screen <- 0
  optRES$print2screen <- 0
  optUNR$CalcSE <- 0
  optRES$CalcSE <- 0

  mBS  <- FCVARestn(x, k, r, optRES)
  mUNR <- FCVARestn(x, k, r, optUNR)

  cat(sprintf('\nHypothesis test to bootstrap:\n'))
  # cat(H)
  H <- HypoTest(mUNR, mBS)

  # How often should the number of iterations be displayed
  show_iters <- 10


  for (j in 1:B) {


    # Display iteration count every 100 Bootstraps
    if(round((j+1)/show_iters) == (j+1)/show_iters) {
      cat(sprintf('iteration: %1.0f\n', j))
    }


    # (1) generate bootstrap DGP under the null
    xBS <- FCVARsimBS(data, mBS, T)
    # append initial values to bootstrap sample
    BSs <- rbind(data, xBS)

    # (2) estimate unrestricted model
    mUNRbs <-  FCVARestn(BSs, k, r, optUNR)

    # (3) estimate restricted model (under the null)
    mRES <-  FCVARestn(BSs, k, r, optRES)

    # (4) calculate test statistic
    LR[j] <- -2*(mRES$like - mUNRbs$like)


  }


  # Return sorted LR stats
  LRbs <- LR[order(LR)]
  # No need to sort if you count the extreme realizations ( but it looks pretty).

  # Calculate Bootstrap P-value (see ETM p.157 eq 4.62)
  H$pvBS <- sum(LRbs > H$LRstat)/B

  # Print output
  cat(sprintf('Bootstrap results:'))
  cat(sprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n',
              H$loglikUNR, H$loglikR))
  cat(sprintf('Test results (df <- %1.0f):\nLR statistic: \t %3.3f\nP-value: \t %1.3f\n',
              H$df, H$LRstat, H$pv))
  cat(sprintf('P-value (BS): \t %1.3f\n', H$pvBS))


  # Return a list of bootstrap test results.
  FCVARboot_out <- list(
    LRbs = LRbs,
    H = H,
    mBS = mBS,
    mUNR = mUNR
  )

  return(FCVARboot_out)
}


################################################################################
# Define function to generate a distribution of a likelihood ratio
#           test statistic for the rank test.
################################################################################
#
# function [LRbs, H, mBS, mUNR] <- FCVARbootRank(x, k, opt, r1, r2, B)
# Written by Michal Popiel and Morten Nielsen (This version 08.06.2015)
#
# DESCRIPTION: This function generates a distribution of a likelihood ratio
#           test statistic for the rank test using a Wild bootstrap,
#			following the method of Cavaliere, Rahbek, and Taylor (2010). It
#           takes the two ranks as inputs to estimate the model under the
#           null and the model under the alternative.
#
# Input <- x  (data - if k>0, actual data is used for initial values)
#         k  (number of lags)
#		  opt(estimation options)
#         r1 (rank under the null)
#         r2 (rank under the alternative)
#         B  (number of bootstrap samples)
#
# Output <- LRbs (B x 1 vector simulated likelihood ratio statistics)
#          pv (approximate p-value for LRstat based on bootstrap
#                                                       distribution)
#          H (a Matlab structure containing LR test results, it is
#               identical to the output from HypoTest, with one addition,
#               namely H$pvBS which is the Bootstrap P-value)
#          mBS  (model estimates under the null)
#          mUNR (model estimates under the alternative)
#
################################################################################


FCVARbootRank <- function(x, k, opt, r1, r2, B) {

  # Calculate length of sample to generate, adjusting for initial values
  T <- nrow(x) - opt$N

  # Use first k+1 observations for initial values
  data <- x[1:(k+1), ]

  LR <- matrix(0, nrow = B, ncol = 1)

  # Turn off output and calculation of standard errors for faster computation
  opt$print2screen <- 0
  opt$print2screen <- 0

  mBS  <- FCVARestn(x, k, r1, opt)
  mUNR <- FCVARestn(x, k, r2, opt)

  # Initialize H (a list containing LR test results, it is
  #               identical to the output from HypoTest, with one addition,
  #               namely H$pvBS which is the Bootstrap P-value)
  H <- list(LRstat = NA,
            pvBS = NA)
  H$LRstat <- -2*(mBS$like - mUNR$like)

  for (j in 1:B) {

    print('j = ')
    print(j)

    # Display iteration count every 100 Bootstraps
    if(round((j+1)/10) == (j+1)/10) {
      cat(sprintf('iteration: %1.0f\n', j))
    }

    # print('made it before mBS')
    # (1) generate bootstrap DGP under the null
    xBS <- FCVARsimBS(data, mBS, T)

    # Faiing on the previous line, with zero rank.
    # print('made it after mBS')

    # append initial values to bootstrap sample
    BSs <- rbind(data, xBS)

    # (2) estimate unrestricted model
    mUNRbs <-  FCVARestn(BSs, k, r2, opt)

    # (3) estimate restricted model (under the null)
    mRES <-  FCVARestn(BSs, k, r1, opt)

    # (4) calculate test statistic
    LR[j] <- -2*(mRES$like - mUNRbs$like)


  }


  # Return sorted LR stats
  LRbs <- LR[order(LR)]

  # Calculate Bootstrap P-value (see ETM p.157 eq 4.62)
  H$pvBS <- sum(LRbs > H$LRstat)/B

  # Print output
  cat(sprintf('Bootstrap results:'))
  cat(sprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n',
              mUNR$like, mBS$like))
  cat(sprintf('Test results:\nLR statistic: \t %3.3f\nP-value (BS): \t %1.3f\n',
              H$LRstat, H$pvBS))


  # Return a list of bootstrap test results.
  FCVARbootRank_out <- list(
    LRbs = LRbs,
    H = H,
    mBS = mBS,
    mUNR = mUNR
  )

  return(FCVARbootRank_out)
}










#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------



################################################################################
# Define function to...
################################################################################
#

#
################################################################################



################################################################################
# End
################################################################################
