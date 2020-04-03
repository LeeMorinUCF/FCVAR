


#' Select Lag Order
#'
#' \code{LagSelect} takes a matrix of variables and performs lag
#' 	selection on it by using the likelihood ratio test. Output and test
#' 	results are printed to the screen.
#'
#' @param x A matrix of variables to be included in the system.
#' @param kmax The maximum number of lags in the system.
#' @param r The cointegrating rank.
#' @param order The order of serial correlation for white noise tests.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{EstOptions()}.
#' @return NULL
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' LagSelect(x, kmax = 4, r = 3, order = 12, opt)
#' @family FCVAR specification functions
#' @seealso \code{EstOptions} to set default estimation options.
#' \code{FCVARestn} is called repeatedly within this function
#' for each candidate lag order.
#' @export
#'
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
  opt <- updateRestrictions(opt, p, r)


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


#' Test for Cointegrating Rank
#'
#' \code{RankTests} performs a sequence of  likelihood ratio tests
# 	for cointegrating rank.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{EstOptions()}.
#' @return A list object \code{rankTestStats} containing the results
#' from cointegrating rank tests, containing the following \code{(p+1)} vectors
#' with \code{i}th element corresponding to \code{rank = i-1}:
#' including the following parameters:
#' \describe{
#'   \item{\code{dHat}}{Estimates of \code{d}.}
#'   \item{\code{bHat}}{Estimates of \code{b}.}
#'   \item{\code{LogL}}{Maximized log-likelihood.}
#'   \item{\code{LRstat}}{LR trace statistic for testing rank \code{r} against rank \code{p}.}
#'   \item{\code{pv}}{The p-value of LR trace test, or "999" if p-value is not available.}
#' }
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' RankTests(x, k = 2, opt)
#' @family FCVAR specification functions
#' @seealso \code{EstOptions} to set default estimation options.
#' \code{FCVARestn} is called repeatedly within this function
#' for each candidate cointegrating rank.
#' @export
#'
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


