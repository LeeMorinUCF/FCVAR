
#' Draw Bootstrap Samples from the FCVAR Model
#'
#' \code{FCVARsimBS} simulates the FCVAR model as specified by
#' input \code{model} and starting values specified by \code{data}.
#' It creates a bootstrap sample by augmenting each iteration
#' with a bootstrap error. The errors are sampled from the
#' residuals specified under the \code{model} input and have a
#' positive or negative sign with equal probability (the Rademacher distribution).
#' @param data A \eqn{T \times p} matrix of starting values for the simulated realizations.
#' @param model A list of estimation results, just as if estimated from \code{FCVARest}.
#' The parameters in \code{model} can also be set or adjusted by assigning new values.
#' @param NumPeriods The number of time periods in the simulation.
#' @return A \code{NumPeriods} \eqn{\times p} matrix \code{xBS} of simulated bootstrap values.
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' model <- FCVARestn(x,k = 3,r = 1,opt)
#' data <- x[1:10, ]
#' xBS <- FCVARsimBS(data, model, NumPeriods = 100)
#' @family FCVAR auxilliary functions
#' @seealso \code{EstOptions} to set default estimation options.
#' \code{FCVARestn} for the specification of the \code{model}.
#' @export
#'
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


#' Calculate Parameters of the Cointegrated VAR Model
#'
#' \code{GetParams} calculates the parameters of the cointegrated VAR model
#'  from a multivariate sample.
#'  This function uses FWL and reduced rank regression to obtain
#' the estimates of \code{Alpha}, \code{Beta}, \code{Rho}, \code{Pi},
#' \code{Gamma}, and \code{Omega}, using the notation in Johansen and Nielsen (2012).
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param db The orders of fractional integration.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{EstOptions()}.
#' @return A list object \code{estimates} containing the estimates,
#' including the following parameters:
#' \describe{
#'   \item{\code{db}}{The orders of fractional integration (taken directly from the input).}
#'   \item{\code{alphaHat}}{A \eqn{p x r} matrix of adjustment parameters.}
#'   \item{\code{betaHat}}{A \eqn{p x r} matrix of cointegrating vectors.
#'       The \eqn{r x 1} vector \eqn{\betax_t} is the stationary cointegration relations.}
#'   \item{\code{rhoHat}}{A \eqn{p x 1} vector of restricted constatnts.}
#'   \item{\code{piHat}}{A \eqn{p x p} matrix \eqn{\Pi = \alpha \beta'} of long-run levels.}
#'   \item{\code{OmegaHat}}{A \eqn{p x p} covariance matrix of the error terms.}
#'   \item{\code{GammaHat}}{A ( \eqn{p x kp} matrix \code{cbind(GammaHat1,...,GammaHatk)})
#'   of autoregressive coefficients. }
#' }
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' estimates <- GetParams(x, k = 2, r = 1, db = c(1, 1), opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{EstOptions} to set default estimation options.
#' \code{FCVARestn} calls \code{GetParams} to estimate the FCVAR model.
#' @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
#' "Likelihood inference for a fractionally cointegrated
#' vector autoregressive model," Econometrica 80, 2667-2732.
#' @export
#'
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

  # print('Z2hat = ')
  # print(Z2hat)

  if ( (k == 0) || is.infinite(kappa(t(Z2hat) %*% Z2hat)) ) {

    # No lags, so there are no effects of Z2.
    R0 <- Z0hat
    R1 <- Z1hat
    # Includes the case when Z2hat is near zero and
    # t(Z2hat) %*% Z2hat is ill-conditioned.

  } else {

    # print('dim(Z2hat) = ')
    # print(dim(Z2hat))
    # print('summary(Z2hat) = ')
    # print(summary(Z2hat))
    # print('t(Z2hat) %*% Z2hat = ')
    # print(t(Z2hat) %*% Z2hat)
    # print('kappa(t(Z2hat) %*% Z2hat) = ')
    # print(kappa(t(Z2hat) %*% Z2hat))

    # Lags included: Obtain the residuals from FWL regressions.
    # Unless, of course, Z2hat is near zero and
    # t(Z2hat) %*% Z2hat is ill-conditioned.
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
    # betaStar <- NULL # Doesn't work for PiHat below
    betaStar <- NA
    # alphaHat <- NULL # Doesn't work for PiHat below
    alphaHat <- NA
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

    # V <- V2[1:p1, ]
    V <- V2[1:p1, , drop = FALSE]
    # betaStar <- V( 1:p1, p1 : -1 : p1-r+1 )
    # betaStar <- matrix(V[ 1:p1, seq(p1, p1-r+1, by = -1) ],
    #                    nrow = p1, ncol = p1-r-1)
    betaStar <- V[ 1:p1, seq(p1, p1-r+1, by = -1), drop = FALSE]
    # Error when k = 2 and r = 2:
    # data length exceeds size of matrix



    # If either alpha or beta is restricted, then the likelihood is
    #   maximized subject to the constraints. This section of the code
    #   replaces the call to the switching algorithm in the previous
    #   version.
    if (!is.null(opt$R_Alpha) | !is.null(opt$R_Beta) ) {


      # [ betaStar, alphaHat, OmegaHat ]...
      # <- RstrctOptm_Switch(betaStar, S00, S01, S11, T, p, opt)
      switched_mats <- RstrctOptm_Switch(betaStar, S00, S01, S11, T, p, opt)
      betaStar <- switched_mats$betaStar
      alphaHat <- switched_mats$alphaHat
      OmegaHat <- switched_mats$OmegaHat

      betaHat <- betaStar[1:p, 1:r, drop = FALSE]
      PiHat <- alphaHat %*% t(betaHat)



    } else {

      # print('betaStar = ')
      # print(betaStar)
      # print('S11 = ')
      # print(S11)

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
    } else {
      rhoHat <- NULL
    }


  } else {

    # (r = p) and do no need reduced rank regression
    V <- S01 %*% solve(S11)
    betaHat <- t(V[, 1:p, drop = FALSE])
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

    # print('V = ')
    # print(V)
    # print('betaHat = ')
    # print(betaHat)
    # print('rhoHat = ')
    # print(rhoHat)

    # Check conformability:
    betaStar <- rbind(betaHat, rhoHat)

  }

  # print('betaHat = ')
  # print(betaHat)

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
      } else {
        xiHat <- t( solve(t(Z3) %*% Z3) %*% t(Z3) %*% (Z0 - Z2 %*% t(GammaHat)) )
      }

    } else {
      if(r>0) {
        xiHat <- t( solve(t(Z3) %*% Z3) %*% t(Z3) %*% (Z0 - Z1 %*% t(PiStar)) )
      } else {
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


