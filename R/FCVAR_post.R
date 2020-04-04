


#' Multivariate White Noise Tests
#'
#' \code{mv_wntest} performs multivariate tests for white noise.
#' It performs both the Ljung-Box Q-test and the LM-test on individual series
#' for a sequence of lag lengths.
#'
#' @param x A matrix of variables to be included in the system,
#' typically model residuals.
#' @param maxlag The number of lags for serial correlation tests.
#' @param printResults An indicator to print results to screen.
#' @return A list object \code{mv_wntest_out} containing the test results,
#' including the following parameters:
#' \describe{
#'   \item{\code{Q}}{A 1xp vector of Q statistics for individual series.}
#'   \item{\code{pvQ}}{A 1xp vector of P-values for Q-test on individual series.}
#'   \item{\code{LM}}{A 1xp vector of LM statistics for individual series.}
#'   \item{\code{pvLM}}{A 1xp vector of P-values for LM-test on individual series.}
#'   \item{\code{mvQ}}{A multivariate Q statistic.}
#'   \item{\code{pvMVQ}}{A p-value for multivariate Q-statistic using \code{p^2*maxlag}
#'   degrees of freedom.}
#' }
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' results <- FCVARestn(x,k = 3,r = 1, opt)
#' mv_wntest(x = results$Residuals, maxlag = 12, printResults = 1)
#'
#' mv_wntest(x = rnorm(100), maxlag = 10, printResults = 1)
#'
#' mv_wntest(x = cumsum(rnorm(100)), maxlag = 10, printResults = 1)
#' @family FCVAR postestimation functions
#' @seealso \code{EstOptions} to set default estimation options.
#' \code{FCVARestn} produces the residuals intended for this test.
#' \code{LagSelect} uses this test as part of the lag order selection process.
#' @note
#' The LM test should be consistent for heteroskedastic series, Q-test is not.
#'
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


#' Breusch-Godfrey Lagrange Multiplier Test for Serial Correlation
#'
#' \code{LMtest} performs a Breusch-Godfrey Lagrange Multiplier test
#' for serial correlation.
#'
#' @param x A vector or Tx1 matrix of variables to be tested,
#' typically model residuals.
#' @param q The number of lags for the serial correlation tests.
#' @return A list object \code{LMtest_out} containing the test results,
#' including the following parameters:
#' \describe{
#'   \item{\code{LM}}{The LM statistic for individual series.}
#'   \item{\code{pv}}{The p-value for LM-test on individual series.}
#' }
#' @examples
#' LMtest(x = rnorm(100), q = 10)
#'
#' LMtest(x = cumsum(rnorm(100)), q = 10)
#' @family FCVAR postestimation functions
#' @seealso \code{mv_wntest} calls this function to test residuals
#' from the estimation results of \code{FCVARestn}.
#' An alternative test is the Ljung-Box Q-test in \code{Qtest}.
#' @note
#' The LM test is consistent for heteroskedastic series.
#'
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


#' Ljung-Box Q-test for Serial Correlation
#'
#' \code{Qtest} performs a (multivariate) Ljung-Box Q-test for serial correlation, see
# 	Luetkepohl (2005, New Introduction to Multiple Time Series Analysis, p. 169).
#'
#' @param x A vector or Tx1 matrix of variables to be tested,
#' typically model residuals.
#' @param maxlag The number of lags for the serial correlation tests.
#' @return A list object \code{LMtest_out} containing the test results,
#' including the following parameters:
#' \describe{
#'   \item{\code{Qstat}}{A 1xp vector of Q statistics for individual series.}
#'   \item{\code{pv}}{A 1xp vector of P-values for Q-test on individual series.}
#' }
#' @examples
#' Qtest(x = rnorm(100), maxlag = 10)
#'
#' Qtest(x = cumsum(rnorm(100)), maxlag = 10)
#' @family FCVAR postestimation functions
#' @seealso \code{mv_wntest} calls this function to test residuals
#' from the estimation results of \code{FCVARestn}.
#' An alternative test is the Breusch-Godfrey Lagrange Multiplier Test in \code{LMtest}.
#' @note
#' The LM test in \code{LMtest} is consistent for heteroskedastic series,
#' while the Q-test is not.
#' @references H. Luetkepohl (2005) "New Introduction to Multiple Time Series Analysis," Springer, Berlin.
#'
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


