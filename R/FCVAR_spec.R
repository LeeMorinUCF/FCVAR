


#' Select Lag Order
#'
#' \code{FCVARlagSelect} takes a matrix of variables and performs lag
#' 	selection on it by using the likelihood ratio test. Output and test
#' 	results are printed to the screen.
#'
#' @param x A matrix of variables to be included in the system.
#' @param kmax The maximum number of lags in the system.
#' @param r The cointegrating rank.This is often set equal to \code{p},
#' the number of variables in the system, since it is better to overspecify
#' than underspecify the model.
#' @param order The order of serial correlation for white noise tests.
#' @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return An S3 object of type \code{FCVAR_lags} containing the results
#' from repeated estimation of the FCVAR model with different orders
#' of the autoregressive lag length.
#' Note that row \code{j} of each of the vectors in the \code{FCVAR_lags} object
#' contains the associated results for lag length \code{j+1}.
#' The \code{FCVAR_lags} object includes the following parameters:
#' \describe{
#'   \item{\code{D}}{A (\code{kmax} + 1) x 2 vector of estimates of d and b.}
#'   \item{\code{loglik}}{A (\code{kmax} + 1) x 1 vector of log-likelihood values.}
#'   \item{\code{LRtest}}{A (\code{kmax} + 1) x 1 vector of likelihood ratio test statistics for tests of significance of \eqn{\Gamma_{j+1}}.}
#'   \item{\code{pvLRtest}}{A (\code{kmax} + 1) x 1 vector of P-values for the likelihood ratio tests of significance of \eqn{\Gamma_{j+1}}.}
#'   \item{\code{i_aic}}{The lag corresponding to the minimum value of the Akaike information criteria.}
#'   \item{\code{aic}}{A (\code{kmax} + 1) x 1 vector of values of the Akaike information criterion.}
#'   \item{\code{i_bic}}{The lag corresponding to the minimum value of the Bayesian information criteria.}
#'   \item{\code{bic}}{A (\code{kmax} + 1) x 1 vector of values of the Bayesian information criterion.}
#'   \item{\code{pvMVq}}{A scalar P-value for the Q-test for multivariate residual white noise.}
#'   \item{\code{pvWNQ}}{A (\code{kmax} + 1) x 1 vector of P-values for the Q-tests for univariate residual white noise.}
#'   \item{\code{pvWNLM}}{A (\code{kmax} + 1) x 1 vector of P-values for the LM-tests for univariate residual white noise.}
#'   \item{\code{kmax}}{The maximum number of lags in the system.}
#'   \item{\code{r}}{The cointegrating rank. This is often set equal to \code{p},
#'     the number of variables in the system, since it is better to overspecify
#'     than underspecify the model.}
#'   \item{\code{p}}{The number of variables in the system.}
#'   \item{\code{cap_T}}{The sample size.}
#'   \item{\code{order}}{The order of serial correlation for white noise tests.}
#'   \item{\code{opt}}{An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
#'     generated from \code{FCVARoptions()}.}
#' }
#' @examples
#' \donttest{
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' FCVARlagSelectStats <- FCVARlagSelect(x, kmax = 3, r = 3, order = 12, opt)
#' }
#' @family FCVAR specification functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} is called repeatedly within this function
#' for each candidate lag order.
#' \code{summary.FCVAR_lags} prints a summary of the output of \code{FCVARlagSelect} to screen.
#' @export
#'
FCVARlagSelect <- function(x, kmax, r, order, opt ) {


  # Determine (initial) dimensions of system.
  cap_T <- nrow(x) # Length of sample (before truncation for initial values).
  p <- ncol(x) # Number of variables.

  # Do not print output for each WN test.
  printWN <- 0

  # Do not print FCVAR estimation for each lag in the loop.
  print2screen <- opt$print2screen
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

    message(sprintf('Estimating for k = %d and r = %d.', k, r))

    # ----- Estimation ---------#
    results <- FCVARestn(x, k, r, opt)


    message(sprintf('Finished Estimation for k = %d and r = %d.', k, r))


    # ----- Record relevant output ---------%
    loglik[k+1] <- results$like
    D[k+1, ]    <- results$coeffs$db
    aic[k+1]    <- -2*loglik[k+1] + 2*results$fp
    bic[k+1]    <- -2*loglik[k+1] + results$fp*log(cap_T - opt$N)

    # ----- White noise tests ---------%
    MVWNtest_stats <- MVWNtest(results$Residuals, order, printWN)
    pvWNQ[k+1, ] <- MVWNtest_stats$pvQ
    pvWNLM[k+1, ] <- MVWNtest_stats$pvLM
    pvMVq[k+1, ] <- MVWNtest_stats$pvMVQ


    # ----- LR test of lag <- k vs lag <- k-1 -----%
    if (k > 0) {
      LRtest[k+1]   <- 2*(loglik[k+1] - loglik[k])
      pvLRtest[k+1] <- 1 - stats::pchisq(LRtest[k+1], p^2)
    }



  }


  # Find lag corresponding to min of information criteria
  i_aic <- which.min(aic)
  i_bic <- which.min(bic)


  # Return list of lag selection statistics.
  FCVARlagSelectStats <- list(
    D = D,
    loglik = loglik,
    LRtest = LRtest,
    pvLRtest = pvLRtest,
    i_aic = i_aic,
    aic = aic,
    i_bic = i_bic,
    bic = bic,
    pvMVq = pvMVq,
    pvWNQ = pvWNQ,
    pvWNLM = pvWNLM,
    kmax = kmax,
    r = r,
    p = p,
    cap_T = cap_T,
    order = order,
    opt = opt
  )
  class(FCVARlagSelectStats) <- 'FCVAR_lags'

  # Print output if required, restoring original settings.
  opt$print2screen <- print2screen
  if (opt$print2screen) {
    summary(object = FCVARlagSelectStats)
  }

  return(FCVARlagSelectStats)

}

#' Summarize Statistics from Lag Order Selection
#'
#' \code{summary.FCVAR_lags} prints a summary of the table of statistics from
#' the output of \code{FCVARlagSelect}.
#' \code{FCVARlagSelect} takes a matrix of variables and performs lag
#' 	selection on it by using the likelihood ratio test.
#'
#' @param object An S3 object of type \code{FCVAR_lags} containing the results
#' from repeated estimation of the FCVAR model with different orders
#' of the autoregressive lag length. It is the output of \code{FCVARlagSelect}.
#' @param ... additional arguments affecting the summary produced.
#' @return NULL
#' @examples
#' \donttest{
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' FCVAR_lag_1 <- FCVARlagSelect(x, kmax = 3, r = 3, order = 12, opt)
#' summary(object = FCVAR_lag_1)
#' }
#' @family FCVAR specification functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} is called repeatedly within this function
#' for each candidate lag order.
#' \code{summary.FCVAR_lags} prints a summary of the output of \code{FCVARlagSelect} to screen.
#' @export
#'
summary.FCVAR_lags <- function(object, ...) {



  #--------------------------------------------------------------------------------
  # Print output
  #--------------------------------------------------------------------------------

  # create a variable for output strings
  yesNo <- c('No','Yes') # Ironic order, No?

  cat(sprintf('\n--------------------------------------------------------------------------------\n'))
  cat(sprintf('                        Lag Selection Results \n'))
  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf('Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n',
              object$p, object$cap_T))
  cat(sprintf('Order for WN tests:   %6.0f     Number of observations for estimation:  %6.0f \n',
              object$order, object$cap_T - object$opt$N))
  cat(sprintf('Restricted constant:  %6s     Initial values:                         %6.0f\n',
              yesNo[object$opt$rConstant+1], object$opt$N )   )
  cat(sprintf('Unrestricted constant:%6s     Level parameter:                        %6s\n',
              yesNo[object$opt$unrConstant+1], yesNo[object$opt$levelParam+1] ))
  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf('Parameter Estimates and Information Criteria:\n'))
  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf(' k  r    d    b      LogL     LR    pv    AIC       BIC'))


  cat(sprintf('\n'))
  for (k in seq(object$kmax, 0, by = -1) ) {

    cat(sprintf('%2.0f %2.0f %4.3f %4.3f %7.2f %6.2f %5.3f %8.2f',
                k, object$r, object$D[k+1, 1], object$D[k+1, 2], object$loglik[k+1],
                object$LRtest[k+1], object$pvLRtest[k+1], object$aic[k+1]))
    # For AIC add asterisk if min value
    if(k+1 == object$i_aic) {cat(sprintf('*'))} else {cat(sprintf(' '))}
    # Print BIC information criteria and add asterisk if min value
    cat(sprintf(' %8.2f', object$bic[k+1]))
    if(k+1 == object$i_bic) {cat(sprintf('*'))} else {cat(sprintf(' '))}


    cat(sprintf('\n'))

  }

  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf('Tests for Serial Correlation of Residuals: \n'))
  cat(sprintf('--------------------------------------------------------------------------------\n'))

  cat(sprintf(' k   pmvQ'))
  for (i in 1:object$p) {
    cat(sprintf('  pQ%1.0f   pLM%1.0f', i, i))
  }

  cat(sprintf('\n'))
  for (k in seq(object$kmax, 0, by = -1) ) {

    cat(sprintf('%2.0f ', k))

    # Print multivariate white noise test P-values
    cat(sprintf('  %4.2f', object$pvMVq[k+1, ]))
    # Print the individual series white noise test P-values
    for (i in 1:object$p) {
      cat(sprintf('  %4.2f  %4.2f', object$pvWNQ[k+1,i], object$pvWNLM[k+1,i]))
    }

    cat(sprintf('\n'))

  }

  cat(sprintf('--------------------------------------------------------------------------------\n'))

}


#' Test for Cointegrating Rank
#'
#' \code{FCVARrankTests} performs a sequence of  likelihood ratio tests
#' 	for cointegrating rank.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return An S3 object of type \code{FCVAR_ranks} containing the results
#' from cointegrating rank tests, containing the following \code{(p+1)} vectors
#' with \code{i}th element corresponding to \code{rank = i-1},
#' including the following parameters:
#' \describe{
#'   \item{\code{dHat}}{Estimates of \code{d}.}
#'   \item{\code{bHat}}{Estimates of \code{b}.}
#'   \item{\code{LogL}}{Maximized log-likelihood.}
#'   \item{\code{LRstat}}{LR trace statistic for testing rank \code{r} against rank \code{p}.}
#'   \item{\code{pv}}{The p-value of LR trace test, or "999" if p-value is not available.}
#'   \item{\code{k}}{The number of lags in the system.}
#'   \item{\code{p}}{The number of variables in the system.}
#'   \item{\code{cap_T}}{The sample size.}
#'   \item{\code{opt}}{An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
#'     generated from \code{FCVARoptions()}.}
#' }
#' @examples
#' \donttest{
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' rankTestStats <- FCVARrankTests(x, k = 2, opt)
#' }
#' @family FCVAR specification functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} is called repeatedly within this function
#' for each candidate cointegrating rank.
#' \code{summary.FCVAR_ranks} prints a summary of the output of \code{FCVARrankTests} to screen.
#' @export
#'
FCVARrankTests <- function(x, k, opt) {

  cap_T <- nrow(x) - opt$N
  p <- ncol(x)

  # Store user specified options for printing to screen because it will be
  # turned off while looping over ranks.
  tempPrint2Screen <- opt$print2screen

  # Create output storage to be filled.
  bHat   <- matrix(0, nrow = p+1, ncol = 1)
  dHat   <- matrix(0, nrow = p+1, ncol = 1)
  LogL   <- matrix(0, nrow = p+1, ncol = 1)
  LRstat <- matrix(0, nrow = p+1, ncol = 1)
  pv     <- matrix(NA, nrow = p+1, ncol = 1)

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

    message(sprintf('Estimating for k = %d and r = %d.', k, r))

    results <- FCVARestn(x, k, r, opt)

    message(sprintf('Finished Estimation for k = %d and r = %d.', k, r))

    dHat[r+1] <- results$coeffs$db[1]
    bHat[r+1] <- results$coeffs$db[2]
    LogL[r+1] <- results$like
  }




  # Calculate the LR statistics and P-values
  for (r in 0 : (p-1)) {

    LRstat[r+1] <-  - 2*( LogL[r+1] - LogL[p+1] )

    p_val <- NULL
    # Get P-values, if
    # (1) no deterministic terms, or
    # (2) there is only restricted constant and d=b, or
    # (3) there is only a level parameter and d=b.
    # In all cases, p-values are only calculated
    #   rank of system (iq) is an integer from 1 through 12.
    if (bHat[r+1] > 0 & bHat[r+1] < 2 & (
      (!opt$rConstant & !opt$unrConstant & !opt$levelParam) |
      (opt$rConstant  & !opt$unrConstant & opt$restrictDB) |
      (opt$levelParam & !opt$unrConstant & opt$restrictDB) )  ) {

      # Call function from fracdist package for p-values.
      if (p - r <= 12) {
        p_val <- fracdist::fracdist_values(iq = p - r,
                                           iscon = consT,
                                           bb = bHat[r+1],
                                           stat = LRstat[r+1])
      } else {
        warning(sprintf('P-values not calculated for the rank test with rank %d.\n', r),
                'Tables of simulated quantiles are only avaiable for rank 12 or lower.',
                'P-values can be calculated by simulation with the FCVARbootRank() function.')
      }


    } else {
      warning(sprintf('P-values not calculated for the rank test with rank %d.\n', r),
              'P-values are only calculated if:\n',
              '1. there are no deterministic terms, or\n',
              '2. there is only restricted constant and d = b, or\n',
              '3. there is only a level parameter and d = b.\n')
    }


    # Store P-values, if calculated.
    if(!is.null(p_val)) {
      pv[r+1] <- p_val
    }

  }



  # Return list of rank test results.
  rankTestStats <- list(
    dHat   = dHat,
    bHat   = bHat,
    LogL   = LogL,
    LRstat = LRstat,
    pv     = pv,
    k      = k,
    p      = p,
    cap_T  = cap_T,
    opt    = opt
  )
  class(rankTestStats) <- 'FCVAR_ranks'


  # Restore settings.
  opt$print2screen <- tempPrint2Screen

  # Print the results to screen.
  if (opt$print2screen) {

    summary(object = rankTestStats)

  }

  return(rankTestStats)
}


#' Summarize Results of Tests for Cointegrating Rank
#'
#' \code{summary.FCVAR_ranks} prints the table of statistics from
#' the output of \code{FCVARrankTests}.
#' \code{FCVARrankTests} performs a sequence of  likelihood ratio tests
#' 	for cointegrating rank.
#'
#' @param object An S3 object of type \code{FCVAR_ranks} containing the results
#' from repeated estimation of the FCVAR model with different
#' cointegrating ranks. It is the output of \code{FCVARrankTests}.
#' @param ... additional arguments affecting the summary produced.
#' @return NULL
#' @examples
#' \donttest{
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' rankTestStats <- FCVARrankTests(x, k = 2, opt)
#' summary(object = rankTestStats)
#' }
#' @family FCVAR specification functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} is called repeatedly within this function
#' for each candidate cointegrating rank.
#' \code{summary.FCVAR_ranks} prints a summary of the output of \code{FCVARrankTests} to screen.
#' @export
#'
summary.FCVAR_ranks <- function(object, ...) {


  # create a variable for output strings
  yesNo <- c('No','Yes')

  cat(sprintf('\n--------------------------------------------------------------------------------\n'))
  cat(sprintf('             Likelihood Ratio Tests for Cointegrating Rank                               \n'))
  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf('Dimension of system:  %6.0f     Number of observations in sample:       %6.0f \n',
              object$p, object$cap_T + object$opt$N))
  cat(sprintf('Number of lags:       %6.0f     Number of observations for estimation:  %6.0f \n',
              object$k, object$cap_T))
  cat(sprintf('Restricted constant:  %6s     Initial values:                         %6.0f\n',
              yesNo[object$opt$rConstant+1], object$opt$N ))
  cat(sprintf('Unestricted constant: %6s     Level parameter:                        %6s\n',
              yesNo[object$opt$unrConstant+1], yesNo[object$opt$levelParam+1] ))
  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf('Rank     d      b     Log-likelihood   LR statistic   P-value\n'))
  for (i in 1:object$p) {
    if (!is.na(object$pv[i])) {
      cat(sprintf('%2.0f     %5.3f  %5.3f  %15.3f  %13.3f  %8.3f\n',
                  i-1, object$dHat[i], object$bHat[i], object$LogL[i], object$LRstat[i], object$pv[i]))
    }
    else {
      cat(sprintf('%2.0f     %5.3f  %5.3f  %15.3f  %13.3f      ----\n',
                  i-1, object$dHat[i], object$bHat[i], object$LogL[i], object$LRstat[i]))
    }

  }

  cat(sprintf('%2.0f     %5.3f  %5.3f  %15.3f           ----      ----\n',
              i, object$dHat[i+1], object$bHat[i+1], object$LogL[i+1]))
  cat(sprintf('--------------------------------------------------------------------------------\n'))


}


#' Distribution of LR Test Statistic for the Rank Test
#'
#' \code{FCVARbootRank} generates a distribution of a likelihood ratio
#'  test statistic for the rank test using a wild bootstrap,
#'	following the method of Cavaliere, Rahbek, and Taylor (2010). It
#'  takes the two ranks as inputs to estimate the model under the
#'  null and the model under the alternative.
#' @param x A matrix of variables to be included in the system.
#' If \code{k>0}, actual data is used for initial values.
#' @param k The number of lags in the system.
#' @param r1 The cointegrating rank under the null hypothesis.
#' @param r2 The cointegrating rank under the alternative hypothesis.
#' @param B The number of bootstrap samples.
#' @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A list \code{FCVARbootRank_stats} containing the test results,
#' including the following parameters:
#' \describe{
#'   \item{\code{LRbs}}{A B x 1 vector of simulated likelihood ratio statistics.}
#'   \item{\code{pv}}{An approximate p-value for the LR statistic based on the bootstrap distribution. }
#'   \item{\code{H}}{A list containing LR test results. It is
#'   identical to the output from \code{HypoTest}, with one addition,
#'   namely \code{H$pvBS} which is the bootstrap p-value)}
#'   \item{\code{mBS}}{Model estimates under the null hypothesis. }
#'   \item{\code{mUNR}}{Model estimates under the alternative hypothesis. }
#' }
#' @examples
#' \donttest{
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$plotRoots <- 0
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' set.seed(42)
#' FCVARbootRank_stats <- FCVARbootRank(x, k = 2, opt, r1 = 0, r2 = 1, B = 2)
#' # In practice, set the number of bootstraps so that (B+1)*alpha is an integer,
#' # where alpha is the chosen level of significance.
#' # For example, set B = 999 (but it takes a long time to compute).
#' }
#' @family FCVAR specification functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{HypoTest} for the format of a hypothesis test results.
#' \code{FCVARestn} for the estimates from a rectricted and unrestricted model within a hypothesis test.
#' @export
#' @references Cavaliere, G., A. Rahbek, and A. M. R. Taylor (2010).
#' "Testing for co-integration in vector autoregressions
#' with non-stationary volatility," Journal of Econometrics 158, 7-24.
#'
FCVARbootRank <- function(x, k, opt, r1, r2, B) {

  # Calculate length of sample to generate, adjusting for initial values
  cap_T <- nrow(x) - opt$N

  # Use first k+1 observations for initial values
  data <- x[1:(k+1), ]

  LR <- matrix(0, nrow = B, ncol = 1)

  # Turn off output and calculation of standard errors for faster computation
  print2screen <- opt$print2screen
  opt$print2screen <- 0
  opt$CalcSE <- 0
  opt$plotRoots <- 0

  mBS  <- FCVARestn(x, k, r1, opt)
  mUNR <- FCVARestn(x, k, r2, opt)

  # Initialize H (a list containing LR test results, it is
  #               identical to the output from HypoTest, with one addition,
  #               namely H$pvBS which is the Bootstrap P-value)
  H <- list(LRstat = NA,
            pvBS = NA)
  H$LRstat <- -2*(mBS$like - mUNR$like)

  # How often should the number of iterations be displayed
  show_iters <- 10

  for (j in 1:B) {

    # Display replication count every show_iters Bootstraps
    if(round((j+1)/show_iters) == (j+1)/show_iters) {
      # cat(sprintf('iteration: %1.0f\n', j))
      message(sprintf('Completed bootstrap replication %d of %d.', j, B))
    }

    # (1) generate bootstrap DGP under the null
    xBS <- FCVARsimBS(data, mBS, cap_T)

    # Append initial values to bootstrap sample
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

  # Print output, if required, after restoring settings.
  opt$print2screen <- print2screen
  if (opt$print2screen) {
    cat(sprintf('Bootstrap rank test results:'))
    cat(sprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n',
                mUNR$like, mBS$like))
    cat(sprintf('Test results:\nLR statistic: \t %3.3f\nP-value (BS): \t %1.3f\n',
                H$LRstat, H$pvBS))
  }


  # Return a list of bootstrap test results.
  FCVARbootRank_stats <- list(
    LRbs = LRbs,
    H = H,
    mBS = mBS,
    mUNR = mUNR
  )

  return(FCVARbootRank_stats)
}


