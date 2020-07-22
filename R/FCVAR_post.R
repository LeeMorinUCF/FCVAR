

#' Test of Restrictions on FCVAR Model
#'
#' \code{FCVARhypoTest} performs a likelihood ratio test of the null
#' 	hypothesis: "model is \code{modelR}" against the alternative hypothesis:
#' 	"model is \code{modelUNR}".
#'
#' @param modelUNR A list of estimation results created for the unrestricted model.
#' @param modelR A list of estimation results created for the restricted model.
#' @return A list object \code{LRtest} containing the test results,
#' including the following parameters:
#' \describe{
#'   \item{\code{loglikUNR}}{The log-likelihood for the unrestricted model.}
#'   \item{\code{loglikR}}{The log-likelihood for the restricted model.}
#'   \item{\code{df}}{The degrees of freedom for the test.}
#'   \item{\code{LRstat}}{The likelihood ratio test statistic.}
#'   \item{\code{p_LRtest}}{The p-value for the likelihood ratio test.}
#' }
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' m1 <- FCVARestn(x, k = 2, r = 1, opt)
#' opt1 <- opt
#' opt1$R_psi <- matrix(c(1, 0), nrow = 1, ncol = 2)
#' opt1$r_psi <- 1
#' m1r1 <- FCVARestn(x, k = 2, r = 1, opt1)
#' Hdb <- FCVARhypoTest(modelUNR = m1, modelR = m1r1)
#'
#' opt1 <- opt
#' opt1$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
#' m1r2 <- FCVARestn(x, k = 2, r = 1, opt1)
#' Hbeta1 <- FCVARhypoTest(m1, m1r2)
#'
#' opt1 <- opt
#' opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
#' m1r4 <- FCVARestn(x, k = 2, r = 1, opt1)
#' Halpha2 <- FCVARhypoTest(m1, m1r4)
#' @family FCVAR postestimation functions
#' @seealso The test is calculated using the results of two calls to
#' \code{FCVARestn}, under the restricted and unrestricted models.
#' Use \code{FCVARoptions} to set default estimation options for each model,
#' then set restrictions as needed before \code{FCVARestn}.
#' @export
#'

FCVARhypoTest <- function(modelUNR, modelR) {


  # Calculate the test statistic.
  LR_test <- 2*(modelUNR$like - modelR$like)

  # Calculate the degrees of freedom by taking the difference in free
  # parameters between the unrestricted and restricted model.
  df <- modelUNR$fp - modelR$fp

  # Calculate the P-value for the test.
  p_LRtest <- 1 - stats::pchisq(LR_test, df)

  # Print output.
  cat(sprintf('Likelihood ratio test results:'))
  cat(sprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n',
              modelUNR$like, modelR$like))
  cat(sprintf('Test results (df = %1.0f):\nLR statistic: \t %3.3f\nP-value: \t %1.3f\n',
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


#' Bootstrap Likelihood Ratio Test
#'
#' \code{FCVARboot} generates a distribution of a likelihood ratio
#' test statistic using a Wild bootstrap, following the method of
#' Boswijk, Cavaliere, Rahbek, and Taylor (2016). It takes two sets
#' of options as inputs to estimate the model under the null and the
#' unrestricted model.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param optRES A list object that stores the chosen estimation options
#'   for the restricted model, as generated from \code{FCVARoptions()},
#'   with adjustments as necessary.
#' @param optUNR A list object that stores the chosen estimation options
#'   for the unrestricted model.
#' @param B The number of bootstrap samples.
#' @return A list object \code{FCVARboot_stats} containing the estimation results,
#' including the following parameters:
#' \describe{
#'   \item{\code{LRbs}}{A \eqn{B x 1} vector of simulated likelihood ratio statistics}
#'   \item{\code{pv}}{An approximate p-value for the likelihood ratio statistic
#'    based on the bootstrap distribution.}
#'   \item{\code{H}}{A list containing the likelihood ratio test results.
#'   It is identical to the output from \code{FCVARhypoTest}, with one addition,
#'   namely \code{H$pvBS} which is the bootstrap p-value}
#'   \item{\code{mBS}}{The model estimates under the null hypothesis.}
#'   \item{\code{mUNR}}{The model estimates under the alternative hypothesis.}
#' }
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' opt$plotRoots <- 0
#' optUNR <- opt
#' optRES <- opt
#' optRES$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
#' set.seed(42)
#' \dontrun{FCVARboot_stats <- FCVARboot(x, k = 2, r = 1, optRES, optUNR, B = 999)}
#' FCVARboot_stats <- FCVARboot(x, k = 2, r = 1, optRES, optUNR, B = 5)
#' @family FCVAR postestimation functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} is called to estimate the models under the null and alternative hypotheses.
#' @references Boswijk, Cavaliere, Rahbek, and Taylor (2016)
#' "Inference on co-integration parameters in heteroskedastic
#' vector autoregressions," Journal of Econometrics 192, 64-85.
#' @export
#'
FCVARboot <- function(x, k, r, optRES, optUNR, B) {


  # Calculate length of sample to generate, adjusting for initial values
  cap_T <- nrow(x) - optRES$N

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
  H <- FCVARhypoTest(mUNR, mBS)

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
  cat(sprintf('Bootstrap likelihood ratio test results:'))
  cat(sprintf('\nUnrestricted log-likelihood: %3.3f\nRestricted log-likelihood:   %3.3f\n',
              H$loglikUNR, H$loglikR))
  cat(sprintf('Test results (df <- %1.0f):\nLR statistic: \t %3.3f\nP-value: \t %1.3f\n',
              H$df, H$LRstat, H$pv))
  cat(sprintf('P-value (BS): \t %1.3f\n', H$pvBS))


  # Return a list of bootstrap test results.
  FCVARboot_stats <- list(
    LRbs = LRbs,
    H = H,
    mBS = mBS,
    mUNR = mUNR
  )

  return(FCVARboot_stats)
}


#' Forecasts with the FCVAR Model
#'
#' \code{FCVARforecast} calculates recursive forecasts with the FCVAR model.
#'
#' @param x A matrix of variables to be included in the system.
#' The forecast will be calculated using these values as starting values.
#' @param model A list of estimation results, just as if estimated from \code{FCVARest}.
#' The parameters in \code{model} can also be set or adjusted by assigning new values.
#' @param NumPeriods The number of time periods in the simulation.
#' @return A \code{NumPeriods} \eqn{\times p} matrix \code{xf} of forecasted values.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' opt1 <- opt
#' opt1$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
#' m1r4 <- FCVARestn(x, k = 2, r = 1, opt1)
#' xf <- FCVARforecast(x, m1r4, NumPeriods = 12)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} for the specification of the \code{model}.
#' \code{FCVARforecast} calls \code{FracDiff} and \code{Lbk} to calculate the forecast.
#' @export
#'
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

  xf <- x
  for (i in 1:NumPeriods) {


    # Append x with zeros to simplify calculations.
    # x <- rbind(x, matrix(0, nrow = 1, ncol = p))
    # x <- rbind(x, rep(0, p))
    xf <- rbind(xf, rep(0, p))
    cap_T <- nrow(xf)

    # Adjust by level parameter if present.
    if(opt$levelParam) {
      y <- xf - matrix(1, nrow = cap_T, ncol = 1) %*% cf$muHat
    } else {
      y <- xf
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
        z <- z + FracDiff( Lbk(matrix(1, nrow = cap_T, ncol = 1), b, 1), d - b ) %*%
          cf$rhoHat %*% t(cf$alphaHat)
      }

    }


    # Add unrestricted constant if present.
    if(opt$unrConstant) {
      z <- z + matrix(1, nrow = cap_T, ncol = 1) %*% t(cf$xiHat)
    }


    # Add lags if present.
    if(!is.null(cf$GammaHat)) {
      # k <- size(cf$GammaHat,2)/p
      k <- ncol(cf$GammaHat)/p
      z <- z +  FracDiff(  Lbk( y , b, k)  , d) %*% t(cf$GammaHat)
    }


    # Adjust by level parameter if present.
    if(opt$levelParam) {
      z <- z + matrix(1, nrow = cap_T, ncol = 1) %*% cf$muHat
    }


    # Append forecast to x matrix.
    xf <- rbind(xf[1:(cap_T-1),], z[cap_T, ])

  }


  #--------------------------------------------------------------------------------
  # Return forecasts.
  #--------------------------------------------------------------------------------

  # Trim off original data to return forecasts only.
  xf <- xf[(nrow(x)+1):nrow(xf), ]

  return(xf)
}



#' Roots of the Characteristic Polynomial
#'
#' \code{GetCharPolyRoots} calculates the roots of the
#' characteristic polynomial and plots them with the unit circle
#' transformed for the fractional model, see Johansen (2008).
#' \code{print.GetCharPolyRoots} prints the output of
#' \code{GetCharPolyRoots} to screen.
#'
#' @param coeffs A list of coefficients for the FCVAR model.
#' An element of the list of estimation \code{results} output from \code{FCVARestn}.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param p The number of variables in the system.
#' @return A complex vector \code{cPolyRoots} with the roots of the characteristic polynomial.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' cPolyRoots <- GetCharPolyRoots(results$coeffs, opt, k = 2, r = 1, p = 3)
#' @family FCVAR postestimation functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} to estimate the model for which to calculate the roots
#' of the characteristic polynomial.
#' \code{print.GetCharPolyRoots} prints the output of
#' \code{GetCharPolyRoots} to screen.
#' @note The roots are calculated from the companion form of the VAR,
#' where the roots are given as the inverse eigenvalues of the
#' coefficient matrix.
#' @references Johansen, S. (2008). "A representation theory for a class of
#' vector autoregressive models for fractional processes,"
#' Econometric Theory 24, 651-676.
#' @export
#'
GetCharPolyRoots <- function(coeffs, opt, k, r, p) {


  b <- coeffs$db[2]

  # First construct the coefficient matrix for the companion form of the VAR.
  PiStar <- diag(p)
  if (r > 0) {
    PiStar <- PiStar + coeffs$alphaHat %*% t(coeffs$betaHat)
  }

  # print('coeffs$GammaHat = ')
  # print(coeffs$GammaHat)
  # print('i = ')
  # print(i)
  # print('p = ')
  # print(p)
  # print('2:k = ')
  # print(2:k)


  if (k > 0) {
    Gamma1 <- coeffs$GammaHat[ , 1 : p]
    PiStar <- PiStar + Gamma1
    if (k > 1) {
      for (i in 2:k) {

        Gammai <- coeffs$GammaHat[ , seq(((i-1)*p + 1), i*p)]
        GammaiMinus1 <- coeffs$GammaHat[ , seq(((i-2)*p + 1), (i-1)*p)]

        PiStar <- cbind(PiStar, (Gammai - GammaiMinus1))

      }
    }

    Gammak <- coeffs$GammaHat[ , seq(((k-1)*p + 1), k*p)]
    PiStar <- cbind(PiStar, ( - Gammak ))
  }

  # print('PiStar = ')
  # print(PiStar)
  # print('p = ')
  # print(p)
  # print('k = ')
  # print(k)

  # Pad with an identity for the transition of the lagged variables.
  if (k > 0) {
    PiStar <- rbind(PiStar,
                    cbind(diag(p*k),
                          matrix(0, nrow = p*k, ncol = p )))
  }


  # The roots are then the inverse eigenvalues of the matrix PiStar.
  cPolyRoots <- 1 / eigen(PiStar)$values
  cPolyRoots <- cPolyRoots[order(-Mod(cPolyRoots))]

  # Generate graph depending on the indicator plotRoots.
  if (opt$plotRoots) {

    plot.GetCharPolyRoots(cPolyRoots, b, file = NULL, file_ext = NULL)

  }


  return(cPolyRoots)
}


#' Print Roots of the Characteristic Polynomial
#'
#' \code{print.GetCharPolyRoots} prints the output of
#' \code{GetCharPolyRoots} to screen.
#' \code{GetCharPolyRoots} calculates the roots of the
#' characteristic polynomial and plots them with the unit circle
#' transformed for the fractional model, see Johansen (2008).
#'
#' @param cPolyRoots A vector of the roots of the characteristic polynomial.
#' An element of the list of estimation \code{results} output from \code{FCVARestn}.
#' @return NULL
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' cPolyRoots <- GetCharPolyRoots(results$coeffs, opt, k = 2, r = 1, p = 3)
#' \dontrun{print.GetCharPolyRoots(cPolyRoots)}
#' \dontrun{plot.GetCharPolyRoots(cPolyRoots, b = results$coeffs$db[2])}
#' @family FCVAR postestimation functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} to estimate the model for which to calculate the roots
#' of the characteristic polynomial.
#' \code{print.GetCharPolyRoots} prints the output of
#' \code{GetCharPolyRoots} to screen.
#' @note The roots are calculated from the companion form of the VAR,
#' where the roots are given as the inverse eigenvalues of the
#' coefficient matrix.
#' @references Johansen, S. (2008). "A representation theory for a class of
#' vector autoregressive models for fractional processes,"
#' Econometric Theory 24, 651-676.
#' @export
#'
print.GetCharPolyRoots <- function(cPolyRoots) {

  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf(  '    Roots of the characteristic polynomial                                                           \n'))
  cat(sprintf('--------------------------------------------------------------------------------\n'))
  cat(sprintf(  '    Number     Real part    Imaginary part       Modulus                                             \n'))
  cat(sprintf('--------------------------------------------------------------------------------\n'))
  for (j in 1:length(cPolyRoots)) {
    cat(sprintf( '      %2.0f       %8.3f       %8.3f         %8.3f                                        \n',
                 j, Re(cPolyRoots[j]), Im(cPolyRoots[j]), Mod(cPolyRoots[j]) ))
  }

  cat(sprintf('--------------------------------------------------------------------------------\n'))

}

#' Plot Roots of the Characteristic Polynomial
#'
#' \code{plot.GetCharPolyRoots} plots the output of
#' \code{GetCharPolyRoots} to screen or to a file.
#' \code{GetCharPolyRoots} calculates the roots of the
#' characteristic polynomial and plots them with the unit circle
#' transformed for the fractional model, see Johansen (2008).
#'
#' @param cPolyRoots A vector of the roots of the characteristic polynomial.
#' An element of the list of estimation \code{results} output from \code{FCVARestn}.
#' @param b The order of fractional integration.
#' @param file A string path and file name in which to plot a figure.
#' Renders to the plot window if running in RStudio if NULL.
#' @param file_ext A string file extension to indicate the graphics format.
#' Either png or pdf are currently supported.
#' @param xlim The coordinates for the horizontal limits of the axes, passed to \code{plot},
#' otherwise set to double the maximum magnitude of points on the unit circle
#' or the transformed unit circle.
#' @param ylim The coordinates for the vertical limits of the axes, passed to \code{plot},
#' otherwise set to double the maximum magnitude of points on the unit circle
#' or the transformed unit circle.
#' @param main The main title of the plot, passed to \code{plot}.
#' If \code{main == 'default'}, a generic title is used.
#' @return NULL
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' cPolyRoots <- GetCharPolyRoots(results$coeffs, opt, k = 2, r = 1, p = 3)
#' \dontrun{print.GetCharPolyRoots(cPolyRoots)}
#' \dontrun{plot.GetCharPolyRoots(cPolyRoots, b = results$coeffs$db[2])}
#' @family FCVAR postestimation functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} to estimate the model for which to calculate the roots
#' of the characteristic polynomial.
#' \code{print.GetCharPolyRoots} prints the output of
#' \code{GetCharPolyRoots} to screen.
#' @note The roots are calculated from the companion form of the VAR,
#' where the roots are given as the inverse eigenvalues of the
#' coefficient matrix.
#' @references Johansen, S. (2008). "A representation theory for a class of
#' vector autoregressive models for fractional processes,"
#' Econometric Theory 24, 651-676.
#' @export
#'
plot.GetCharPolyRoots <- function(cPolyRoots, b,
                                  file = NULL, file_ext = NULL,
                                  xlim = NULL, ylim = NULL, main = NULL) {

  # print('cPolyRoots = ')
  # print(cPolyRoots)
  # print('b = ')
  # print(b)

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

  # Determine axes based on largest roots, if not specified.
  if (is.null(xlim) & is.null(ylim)) {
    maxXYaxis <- max( c(transformedUnitCircleX, unitCircleX,
                        transformedUnitCircleY, unitCircleY) )
    minXYaxis <- min( c(transformedUnitCircleX, unitCircleX,
                        transformedUnitCircleY, unitCircleY) )
    maxXYaxis <- max( maxXYaxis, -minXYaxis )

    xlim <- 2*c(-maxXYaxis, maxXYaxis)
    ylim <- 2*c(-maxXYaxis, maxXYaxis)
  }

  # Open file, if required.
  if (!is.null(file)) {
    if(file_ext == 'pdf') {
      grDevices::pdf(file)
    } else if(file_ext == 'png') {
      grDevices::png(file)
    } else {
      stop('Graphics format not supported. Try pdf or png format.')
    }

  }

  if (!is.null(main)) {
    if (main == 'default') {
      main <- c('Roots of the characteristic polynomial',
                'with the image of the unit circle')
    }
  }


  graphics::plot(transformedUnitCircleX,
       transformedUnitCircleY,
       main = main,
       xlab = 'Real Part of Root',
       ylab = 'Imaginary Part of Root',
       xlim = xlim,
       ylim = ylim,
       type = 'l',
       lwd = 3,
       col = 'red')
  graphics::lines(unitCircleX, unitCircleY, lwd = 3, col = 'black')
  graphics::points(Re(cPolyRoots), Im(cPolyRoots),
         pch = 16, col = 'blue')

  # Close graphics device, if any.
  if (!is.null(file)) {
    grDevices::dev.off()
  }

}


#' Multivariate White Noise Tests
#'
#' \code{MVWNtest} performs multivariate tests for white noise.
#' It performs both the Ljung-Box Q-test and the LM-test on individual series
#' for a sequence of lag lengths.
#' \code{print.MVWNtest} prints these statistics to screen.
#'
#' @param x A matrix of variables to be included in the system,
#' typically model residuals.
#' @param maxlag The number of lags for serial correlation tests.
#' @param printResults An indicator to print results to screen.
#' @return A list object \code{MVWNtest_stats} containing the test results,
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
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
#'
#' set.seed(27)
#' WN <- stats::rnorm(100)
#' RW <- cumsum(stats::rnorm(100))
#' MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
#' MVWNtest_stats <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)
#' @family FCVAR postestimation functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} produces the residuals intended for this test.
#' \code{LagSelect} uses this test as part of the lag order selection process.
#' \code{print.MVWNtest} prints the output of \code{MVWNtest} to screen.
#' @note
#' The LM test should be consistent for heteroskedastic series, Q-test is not.
#' @export
#'
MVWNtest <- function(x, maxlag, printResults) {

  cap_T <- nrow(x)
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


  # Output a list of results.
  MVWNtest_stats <- list(
    Q = Q,
    pvQ = pvQ,
    LM = LM,
    pvLM = pvLM,
    mvQ = mvQ,
    pvMVQ = pvMVQ
  )

  # Print output
  if (printResults) {

    print.MVWNtest(stats = MVWNtest_stats, maxlag, p)

  }


  return(MVWNtest_stats)
}


#' Print Statistics for Multivariate White Noise Tests
#'
#' \code{print.MVWNtest} prints the statistics from \code{MVWNtest} to screen.
#' \code{MVWNtest} performs multivariate tests for white noise.
#' It performs both the Ljung-Box Q-test and the LM-test on individual series
#' for a sequence of lag lengths.
#'
#' @param stats A list object \code{MVWNtest_stats} containing the results
#' from multivariate tests for white noise
#' It is the output of \code{MVWNtest}.
#' @param maxlag The number of lags for serial correlation tests.
#' @param p The number of variables in the system.
#' @return NULL
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
#' \dontrun{print.MVWNtest(stats = MVWNtest_stats, maxlag = 12, p = 3)}
#'
#' set.seed(27)
#' WN <- stats::rnorm(100)
#' RW <- cumsum(stats::rnorm(100))
#' MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
#' MVWNtest_stats <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)
#' \dontrun{print.MVWNtest(stats = MVWNtest_stats, maxlag = 10, p = 2)}
#' @family FCVAR postestimation functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} produces the residuals intended for this test.
#' \code{LagSelect} uses this test as part of the lag order selection process.
#' \code{print.MVWNtest} prints the output of \code{MVWNtest} to screen.
#' @note
#' The LM test should be consistent for heteroskedastic series, Q-test is not.
#' @export
#'
print.MVWNtest <- function(stats, maxlag, p) {


  cat(sprintf('\n       White Noise Test Results (lag = %g)\n', maxlag))
  cat(sprintf('---------------------------------------------\n'))
  cat(sprintf('Variable |       Q  P-val |      LM  P-val  |\n'))
  cat(sprintf('---------------------------------------------\n'))
  cat(sprintf('Multivar | %7.3f  %4.3f |     ----  ----  |\n', stats$mvQ, stats$pvMVQ))
  for (i in 1:p) {
    cat(sprintf('Var%g     | %7.3f  %4.3f | %7.3f  %4.3f  |\n',
                i, stats$Q[i], stats$pvQ[i], stats$LM[i], stats$pvLM[i] ))
  }

  cat(sprintf('---------------------------------------------\n'))

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
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
#' LMtest(x = matrix(results$Residuals[, 1]), q = 12)
#' LMtest(x = results$Residuals[,2, drop = FALSE], q = 12)
#'
#' set.seed(27)
#' WN <- stats::rnorm(100)
#' RW <- cumsum(stats::rnorm(100))
#' LMtest(x = matrix(WN), q = 10)
#' LMtest(x = matrix(RW), q = 10)
#' MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
#' MVWNtest_stats <- MVWNtest(x = MVWN_x, maxlag = 10, printResults = 1)
#' @family FCVAR postestimation functions
#' @seealso \code{MVWNtest} calls this function to test residuals
#' from the estimation results of \code{FCVARestn}.
#' An alternative test is the Ljung-Box Q-test in \code{Qtest}.
#' @note
#' The LM test is consistent for heteroskedastic series.
#' @export
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
  cap_T <- nrow(x)
  x <- x - mean(x)

  # cat(sprintf('x is %d by %d', nrow(x), ncol(x)))

  # print(x)

  y <- x[seq(q + 1, cap_T), ,  drop = FALSE]
  z <- x[seq(1, cap_T - q), ,  drop = FALSE]

  for (i in 1:(q - 1)) {
    z <- cbind(x[seq(i + 1, cap_T - q + i), ,  drop = FALSE], z)
  }



  # print(head(y))
  # print(tail(y))
  # cat(sprintf('y is %d by %d', nrow(y), ncol(y)))


  e <- y
  # print(head(e))
  # print(tail(e))


  # cat(sprintf('e is %d by %d', nrow(e), ncol(e)))
  # cat(sprintf('z is %d by %d', nrow(z), ncol(z)))
  # kron_test <- kronecker(matrix(1, 1, q), e)
  # cat(sprintf('kron_test is %d by %d', nrow(kron_test), ncol(kron_test)))


  # s <- z[,1:q] * repmat(e,1,q)
  # Translate this properly:
  s <- z[,1:q,  drop = FALSE] * kronecker(matrix(1, 1, q), e)

  # cat(sprintf('s is %d by %d', nrow(s), ncol(s)))

  sbar <- t(colMeans(s))

  # print('sbar = ')
  # print(sbar)

  kron_sbar <- kronecker(matrix(1, nrow(s)), sbar)

  # print(head(kron_sbar))
  # print(tail(kron_sbar))

  # The next line bsxfun(@FUNC, A, B) applies the element-by-element binary
  # operation FUNC to arrays A and B, with singleton expansion enabled.
  # Need to translate this to R:
  # s <- mapply(minus, s, sbar)
  # Never mind that fancy pants stuff for something that
  # is not even calculated in a loop.
  s <- s - kron_sbar

  S <- t(s) %*% s/cap_T
  # LMstat <- T*sbar %*% S^(-1) %*% t(sbar)
  LMstat <- cap_T*sbar %*% solve(S) %*% t(sbar)
  pv <- 1 - stats::pchisq(LMstat, q)


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
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' MVWNtest_stats <- MVWNtest(x = results$Residuals, maxlag = 12, printResults = 1)
#' Qtest(x = results$Residuals, maxlag = 12)
#' Qtest(x = matrix(results$Residuals[, 1]), maxlag = 12)
#' Qtest(x = results$Residuals[,2, drop = FALSE], maxlag = 12)
#'
#' set.seed(27)
#' WN <- stats::rnorm(100)
#' RW <- cumsum(stats::rnorm(100))
#' MVWN_x <- as.matrix(data.frame(WN = WN, RW = RW))
#' Qtest(x = MVWN_x, maxlag = 10)
#' Qtest(x = matrix(WN), maxlag = 10)
#' Qtest(x = matrix(RW), maxlag = 10)
#' @family FCVAR postestimation functions
#' @seealso \code{MVWNtest} calls this function to test residuals
#' from the estimation results of \code{FCVARestn}.
#' An alternative test is the Breusch-Godfrey Lagrange Multiplier Test in \code{LMtest}.
#' @note
#' The LM test in \code{LMtest} is consistent for heteroskedastic series,
#' while the Q-test is not.
#' @references H. Luetkepohl (2005) "New Introduction to Multiple Time Series Analysis," Springer, Berlin.
#' @export
#'
Qtest <- function(x, maxlag) {

  # print('class(x) = ')
  # print(class(x))

  cap_T <- nrow(x)
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
  for (t in 1:cap_T) {
    # C0 <- C0 + x[t, ] %*% t(x[t, ])
    C0 <- C0 + t(x[t, , drop = FALSE]) %*% x[t, , drop = FALSE]
  }
  C0 <- C0/cap_T

  C <- array(rep(0, p*p*maxlag), dim = c(p,p,maxlag))
  for (i in 1:maxlag) {
    for (t in (i + 1):cap_T) {
      C[ , ,i] <- C[ , ,i] + t(x[t, , drop = FALSE]) %*% x[t-i, , drop = FALSE]
    }

    C[ , ,i] <- C[ , ,i]/(cap_T - i) # Note division by (T-i) instead of T.
  }


  # (Multivariate) Q statistic
  Qstat <- 0
  for (j in 1:maxlag) {
    # Qstat <- Qstat+trace(C(:,:,j)'*inv(C0)*C(:,:,j)*inv(C0)) / (T-j) %'
    # The following line is a more efficient calculation than the previous
    # Qstat <- Qstat+trace( (C(:,:,j)'/C0)*(C(:,:,j)/C0) ) / (T-j) %' #'
    # Need function for trace:
    Qstat <- Qstat + sum(diag( (t(C[ , ,j]) %*% solve(C0)) %*% (C[ , ,j] %*% solve(C0)) )) / (cap_T - j)
  }


  Qstat <- Qstat*cap_T*(cap_T + 2)
  pv <- 1 - stats::pchisq(Qstat, p*p*maxlag) # P-value is calculated with p^2*maxlag df.


  # Output a list of results.
  Qtest_out <- list(
    Qstat = Qstat,
    pv = pv
  )

  return(Qtest_out)
}
