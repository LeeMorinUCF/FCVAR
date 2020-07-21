
#' Draw Samples from the FCVAR Model
#'
#' \code{FCVARsim} simulates the FCVAR model as specified by
#' input \code{model} and starting values specified by \code{data}.
#' Errors are drawn from a Normal distribution.
#' @param x A \eqn{N x p} matrix of \code{N} starting values for the simulated observations.
#' @param model A list of estimation results, just as if estimated from \code{FCVARest}.
#' The parameters in \code{model} can also be set or adjusted by assigning new values.
#' @param NumPeriods The number of time periods in the simulation.
#' @return A \code{NumPeriods} \eqn{x p} matrix \code{xBS} of simulated observations.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' x_sim <- FCVARsim(x[1:10, ], results, NumPeriods = 100)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} for the specification of the \code{model}.
#' Use \code{FCVARsim} to draw a sample from the FCVAR model.
#' For simulations intended for bootstrapping statistics, use \code{FCVARsimBS}.
#' @export
#'
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
  err <- matrix(stats::rnorm(NumPeriods*p), nrow = NumPeriods, ncol = p)

  #--------------------------------------------------------------------------------
  # Recursively generate simulated data values
  #--------------------------------------------------------------------------------

  # Initialize simulation with starting values.
  x_sim <- x
  for (i in 1:NumPeriods) {


    # Append x_sim with zeros to simplify calculations.
    x_sim <- rbind(x_sim, rep(0, p))
    T <- nrow(x_sim)

    # Adjust by level parameter if present.
    if(opt$levelParam) {
      y <- x_sim - matrix(1, nrow = T, ncol = 1) %*% cf$muHat
    } else {
      y <- x_sim
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

    # Append to x_sim matrix.
    x_sim <- rbind(x_sim[1:(T-1), ], z[T,])

  }

  #--------------------------------------------------------------------------------
  # Return simulated data values, including (EXcluding?) initial values.
  #--------------------------------------------------------------------------------

  x_sim <- x_sim[(nrow(x)+1):nrow(x_sim), ]

  return(x_sim)
}


#' Draw Bootstrap Samples from the FCVAR Model
#'
#' \code{FCVARsimBS} simulates the FCVAR model as specified by
#' input \code{model} and starting values specified by \code{data}.
#' It creates a bootstrap sample by augmenting each iteration
#' with a bootstrap error. The errors are sampled from the
#' residuals specified under the \code{model} input and have a
#' positive or negative sign with equal probability (the Rademacher distribution).
#' @param data A \eqn{T x p} matrix of starting values for the simulated realizations.
#' @param model A list of estimation results, just as if estimated from \code{FCVARest}.
#' The parameters in \code{model} can also be set or adjusted by assigning new values.
#' @param NumPeriods The number of time periods in the simulation.
#' @return A \code{NumPeriods} \eqn{x p} matrix \code{xBS} of simulated bootstrap values.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' xBS <- FCVARsimBS(x[1:10, ], results, NumPeriods = 100)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} for the specification of the \code{model}.
#' Use \code{FCVARsim} to draw a sample from the FCVAR model.
#' For simulations intended for bootstrapping statistics, use \code{FCVARsimBS}.
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
    2*( matrix(stats::rnorm(T), nrow = T, ncol = 1) > 0) %*% matrix(1, nrow = 1, ncol = p)

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
      z <- z + matrix(1, nrow = T, ncol = 1) %*% t(cf$xiHat)
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
#' generated from \code{FCVARoptions()}.
#' @return A list object \code{estimates} containing the estimates,
#' including the following parameters:
#' \describe{
#'   \item{\code{db}}{The orders of fractional integration (taken directly from the input).}
#'   \item{\code{alphaHat}}{A \eqn{p x r} matrix of adjustment parameters.}
#'   \item{\code{betaHat}}{A \eqn{p x r} matrix of cointegrating vectors.
#'       The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.}
#'   \item{\code{rhoHat}}{A \eqn{p x 1} vector of restricted constatnts.}
#'   \item{\code{piHat}}{A \eqn{p x p} matrix \eqn{\Pi = \alpha \beta'} of long-run levels.}
#'   \item{\code{OmegaHat}}{A \eqn{p x p} covariance matrix of the error terms.}
#'   \item{\code{GammaHat}}{A ( \eqn{p x kp} matrix \code{cbind(GammaHat1,...,GammaHatk)})
#'   of autoregressive coefficients. }
#' }
#' @examples
#' opt <- FCVARoptions()
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' estimates <- GetParams(x, k = 2, r = 1, db = c(1, 1), opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
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
      # <- GetRestrictedParams(betaStar, S00, S01, S11, T, p, opt)
      switched_mats <- GetRestrictedParams(betaStar, S00, S01, S11, T, p, opt)
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


#' Grid Search to Maximize Likelihood Function
#'
#' \code{FCVARlikeGrid} performs a grid-search optimization
#' by calculating the likelihood function
#' on a grid of candidate parameter values.
#' This function evaluates the likelihood over a grid of values
#' 	for \code{c(d,b)} (or \code{phi}).
#' 	It can be used when parameter estimates are sensitive to
#' 	starting values to give an approximation of the global max which can
#' 	then be used as the starting value in the numerical optimization in
#' 	\code{FCVARestn}.
#' 	\code{plot.FCVARlikeGrid} plots the likelihood function from \code{FCVARlikeGrid}.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A list object \code{likeGrid_params} containing the optimization results,
#' including the following parameters:
#' \describe{
#'   \item{\code{params}}{A vector \code{params} of \code{d} and \code{b}
#' (and \code{mu} if level parameter is selected)
#' corresponding to a maximum over the grid of \code{c(d,b)} or \code{phi}.}
#'   \item{\code{dbHatStar}}{A vector of \code{d} and \code{b}
#'   corresponding to a maximum over the grid of \code{c(d,b)} or \code{phi}}
#'   \item{\code{muHatStar}}{A vector of the optimal \code{mu} if level parameter is selected. }
#'   \item{\code{Grid2d}}{An indicator for whether or not the optimization
#'   is conducted over a 2-dimensional parameter space,
#'   i.e. if there is no equality restriction on \code{d} and \code{b}.}
#'   \item{\code{dGrid}}{A vector of the grid points in the parameter \code{d}.}
#'   \item{\code{bGrid}}{A vector of the grid points in the parameter \code{b}.}
#'   \item{\code{like}}{The maximum value of the likelihood function over the chosen grid.}
#' }
#'
#' @examples
#' # Restrict equality of fractional parameters.
#' opt <- FCVARoptions()
#' opt$dbStep1D     <- 0.1 # Coarser grid for plotting example.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$restrictDB   <- 1 # impose restriction d=b ? 1 <- yes, 0 <- no.
#' opt$progress     <- 2 # Show progress report on each value of b.
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
#' \dontrun{plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt, main = 'default')}
#'
#' # Linear restriction on fractional parameters.
#' opt <- FCVARoptions()
#' opt$dbStep1D     <- 0.1 # Coarser grid for plotting example.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
#' # Impose linear restriction on d and b:
#' opt$R_psi        <- matrix(c(2, -1), nrow = 1, ncol = 2)
#' opt$r_psi        <- 0.5
#' opt$progress     <- 2 # Show progress report on each value of b.
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
#' \dontrun{plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt, main = 'default')}
#'
#' # Constrained 2-dimensional optimization.
#' # Impose restriction dbMax >= d >= b >= dbMin.
#' \dontrun{
#' opt <- FCVARoptions()
#' opt$dbStep1D     <- 0.1 # Coarser grid for plotting example.
#' opt$dbStep2D     <- 0.2 # Coarser grid for plotting example.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 1 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
#' opt$progress     <- 2 # Show progress report on each value of b.
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
#' }
#'
#' # Unconstrained 2-dimensional optimization.
#' \dontrun{
#' opt <- FCVARoptions()
#' opt$dbStep1D     <- 0.01 # Coarser grid for plotting example.
#' opt$dbStep2D     <- 0.2 # Coarser grid for plotting example.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$restrictDB   <- 0 # impose restriction d=b ? 1 <- yes, 0 <- no.
#' opt$progress     <- 2 # Show progress report on each value of b.
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
#' }
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{plot.FCVARlikeGrid} plots the likelihood function from \code{FCVARlikeGrid}.
#' @note If \code{opt$LocalMax == 0}, \code{FCVARlikeGrid} returns the parameter values
#'       corresponding to the global maximum of the likelihood on the grid.
#'       If \code{opt$LocalMax == 1}, \code{FCVARlikeGrid} returns the parameter values for the
#'       local maximum corresponding to the highest value of \code{b}. This
#'       alleviates the identification problem mentioned in Johansen and
#'       Nielsen (2010, section 2.3).
#' @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2010).
#' "Likelihood inference for a nonstationary fractional
#' autoregressive model," Journal of Econometrics 158, 51-66.
#' @export
#'
FCVARlikeGrid <- function(x, k, r, opt) {



  p <- ncol(x)

  if(opt$progress != 0) {
    if(opt$progress == 1) {
      # This is cute but... not right now:
      # M_status_bar <- waitbar(0,['Model: k= ',num2str(k), ', r= ', num2str(r)])
    }

    # lastTic <- tic()
  }

  #--------------------------------------------------------------------------------
  # INITIALIZATION
  #--------------------------------------------------------------------------------

  # Check if performing a 2-dimensional or 1-dimensional grid search. The
  #   step size will be smaller if searching in only one dimension.
  if(is.null(opt$R_psi)) {
    Grid2d <- 1
    # dbStep <- 0.02
    # dbStep <- 0.2
    dbStep <- opt$dbStep2D # Make this a variable in options.
  } else {
    Grid2d <- 0
    # dbStep <- 0.01
    dbStep <- opt$dbStep1D
  }



  # If equality restrictions are imposed, need to construct a grid over
  #   phi and obtain db <- H*phi + h.
  if(!is.null(opt$R_psi)) {
    # This set of variables is defined for easy translation between
    #   phi (unrestricted parameter) and d,b (restricted parameters).
    R <- opt$R_psi
    s <- opt$r_psi
    H <- pracma::null(R)
    h <- t(R) %*% solve(R %*% t(R)) %*% s
  }


  # GetBounds returns upper and lower bounds in terms of phi when
  # restrictions are imposed and in terms of d and b when they are
  # unrestricted. It also adjusts for limits on d and b imposed by user.
  UB_LB_bounds <- GetBounds(opt)
  dbMax <- UB_LB_bounds$UB
  dbMin <- UB_LB_bounds$LB



  # Set up the grid.
  if(Grid2d) {

    # Search over d as well as b.
    # Set up grid for d parameter
    dGrid <-  seq(dbMin[1], dbMax[1], by = dbStep)
    # Number of grid points
    nD    <- length(dGrid)
    # Set up grid for b parameter
    bGrid <-  seq(dbMin[2], dbMax[2], by = dbStep)
    # Number of grid points
    nB    <- length(bGrid)
    if(opt$constrained) {
      # Since it is possible to have different grid lengths for d and
      # b when the parameters have different max/min values, the
      # following function calculates the total number of iterations
      # by counting the number of instances in which the grid values
      # in b are less than or equal to those in d. This should be
      # moved to the beginning of the loop for faster computation
      # later.
      # totIters <- sum(sum(bsxfun(@le,repmat(bGrid',1, nD),dGrid))) %' ))))
      bdGrid <- expand.grid(dGrid = dGrid, bGrid = bGrid)
      totIters <- sum(bdGrid[, 'bGrid'] <= bdGrid[, 'dGrid'])
    } else {
      # Unconstrained so search through entire grid for d.
      dStart   <- 1
      totIters <- nB*nD
    }


  } else {
    # Only searching over one parameter.
    bGrid <-  seq(dbMin[1], dbMax[1], by = dbStep)
    dGrid <-  NA
    nB    <- length(bGrid)
    nD    <- 1
    dStart   <- 1
    totIters <- nB
  }



  # Create a matrix of NA's to store likelihoods, we use NA's here
  #   because NA entries are not plotted and do not affect the finding
  #   of the maximum.
  like  <- matrix(NA, nrow = nB, ncol = nD)

  # Initialize storage bin and starting values for optimization involving
  #   level parameter.
  if(opt$levelParam) {
    mu  <- array(0, dim = c(p, nB, nD))
    StartVal <- x[1, ,drop = FALSE]
  }


  #--------------------------------------------------------------------------------
  # CALCULATE LIKELIHOOD OVER THE GRID
  #--------------------------------------------------------------------------------



  iterCount <- 0

  for (iB in 1:nB) {

    b <- bGrid[iB]

    if (opt$progress == 2) {
      message(sprintf('Now estimating for iteration %d of %d: b = %f.',
                      iterCount, totIters, b))
    }

    # If d>=b (constrained) then search only over values of d that are
    #   greater than or equal to b. Also, perform this operation only if
    #   searching over both d and b.
    if(opt$constrained & Grid2d) {
      # Start at the index that corresponds to the first value in the
      # grid for d that is >= b
      # dStart <-  find(dGrid>=b, 1)
      # dStart <- dGrid[dGrid >= b][1]
      dStart <- which(dGrid >= b)[1]
    }


    for (iD in dStart:nD) {

      iterCount <- iterCount + 1

      # db is definied in terms of d,b for use by FCVARlikeMU
      # if level parameters are present and for displaying in
      # the output. phi is used by FCVARlike, which can
      # handle both phi or d,b and makes appropriate
      # adjustments inside the function.
      if(is.null(opt$R_psi)) {
        d <- dGrid[iD]
        db <- cbind(d, b)
        phi <- db
      } else {
        phi <- bGrid[iB]
        db <- H*phi + h
      }


      if(opt$levelParam) {
        # Optimize over level parameter for given (d,b).
        # [ muHat, maxLike, ! ] ...
        # <- fminunc(@( params ) -FCVARlikeMu(x, db, params, k, r, opt), ...
        #            StartVal, opt$UncFminOptions )

        # print('Made it here!')

        min_out <- stats::optim(StartVal,
                         function(params) {-FCVARlikeMu(params, x, db, k, r, opt)})

        # print('Didnt make it here!')

        muHat <- min_out$par
        maxLike <- min_out$value

        # Store the results.
        like[iB,iD] <- -maxLike
        mu[, iB,iD] <- muHat

      } else {
        # Only called if no level parameters. phi contains
        # either one or two parameters, depending on
        # whether or not restrictions have been imposed.
        like[iB,iD] <- FCVARlike(phi, x, k, r, opt)
      }


      if (opt$progress != 0) {
        # if (toc(lastTic) > opt$updateTime | iterCount == totIters) {
        if ( iterCount == totIters) {
          if (opt$progress == 1) {
            SimNotes <- sprintf('Model: k=%g, r=%g\nb=%4.2f, d=%4.2f, like=%8.2f',
                                k, r, db[2],db[1], like[iB,iD] )
            # waitbar(iterCount/totIters,M_status_bar, SimNotes)
          } else {
            message(sprintf('Progress : %5.1f%%, b=%4.2f, d=%4.2f, like=%g.',
                            (iterCount/totIters)*100, db[2], db[1], like[iB,iD] ))
          }

          # lastTic <- tic()
        }

      }


    }


  }

  # # Save the workspace at this point.
  # save.image(file = 'FCVARlikeGridData1.RData')
  #
  # # Testing: Analyze the like grid
  # nrow(like)
  # ncol(like)
  # summary(like)
  # max(like)
  # min(like)
  #
  # # Stack into a data frame.
  # like_df <- data.frame(expand.grid(b = bGrid,
  #                                   d = dGrid),
  #                       like = NA)
  # head(like_df)
  # tail(like_df)
  #
  # for (col_num in 1:ncol(like)) {
  #   like_df[seq((col_num - 1)*length(bGrid) + 1,
  #             col_num*length(bGrid)), 'like'] <- like[, col_num]
  # }
  #
  # head(like_df[order(-like_df[ , 'like']), ])
  # tail(like_df[order(-like_df[ , 'like']), ])
  #
  # # Test some values.
  # like[bGrid == 0.49, dGrid == 0.01]
  # like[bGrid > 0.469 & bGrid < 0.471, dGrid == 0.01]

  # Restricted version:
  # like_df <- data.frame(cbind(bGrid, like))
  # colnames(like_df) <- c('phi', 'like')
  # like_df[, 'b'] <- NA
  # like_df[, 'b'] <- H[1]*like_df[, 'phi'] + h[1]
  # head(like_df[order(-like_df[ , 'like']), ])
  #
  # head(like_df[order(-like_df[ , 'like']), ][like_df[ , 'b'] > 0.65 &
  #                                              like_df[ , 'b'] < 0.68, ])


  # db <- H*phi + h

  #--------------------------------------------------------------------------------
  # FIND THE MAX OVER THE GRID
  #--------------------------------------------------------------------------------


  if(opt$LocalMax) {

    # Local max.
    if(Grid2d) {
      # [!,smax,!,!] <- extrema2(like,1)
      # [indexB,indexD] <- ind2sub(size(like),smax)
      indexB_and_D <- which(like == max(like), arr.ind = TRUE)
      indexB <- indexB_and_D[, 1, drop = FALSE]
      indexD <- indexB_and_D[, 2, drop = FALSE]

    } else {
      # [!,indexB,!,!] <- extrema(like)
      indexB <- which(like == max(like))
      indexD <- 1
    }


    # If there is no local max, return global max.
    if(is.null(indexB) | is.null(indexD)) {
      # [ indexB, indexD ] <- find(like == max(max(like)))

      indexB_and_D <- which(like == max(like), arr.ind = TRUE)
      indexB <- indexB_and_D[, 1, drop = FALSE]
      indexD <- indexB_and_D[, 2, drop = FALSE]
    }


  } else {
    # Global max.
    # [ indexB, indexD ] <- find(like == max(max(like)))

    indexB_and_D <- which(like == max(like), arr.ind = TRUE)
    indexB <- indexB_and_D[, 1, drop = FALSE]
    indexD <- indexB_and_D[, 2, drop = FALSE]
  }


  if(length(indexD) > 1 | length(indexB) > 1) {
    # If maximum is not unique, take the index pair corresponding to
    # the highest value of b.
    if(Grid2d) {
      # Sort in ascending order according to indexB.
      # [indexB,indBindx] <- sort(indexB)
      indexB <- indexB[order(indexB)]
      indBindx <- order(indexB)
      indexD <- indexD[indBindx]
    } else {
      # Sort in ascending order.
      indexB <- indexB[order(indexB)]
    }

    indexB <- indexB[length(indexB)]
    indexD <- indexD[length(indexD)]
    # cat(sprintf('\nWarning, grid search did not find a unique maximum of the log-likelihood function.\n')  )
    warning('Grid search did not find a unique maximum of the log-likelihood function.')
  }

  # Translate to d,b.
  if(!is.null(opt$R_psi)) {
    dbHatStar <- t(H %*% bGrid[indexB] + h)
  } else {
    dbHatStar <- cbind(dGrid[indexD], bGrid[indexB])
  }


  #--------------------------------------------------------------------------------
  # STORE THE PARAMETER VALUES
  #--------------------------------------------------------------------------------

  # Output a list of values, including the previous vector \code{params}.
  likeGrid_params <- list(
    params = NA,
    dbHatStar = dbHatStar,
    muHatStar = NA,
    Grid2d = Grid2d,
    dGrid = dGrid,
    bGrid = bGrid,
    like = like
  )

  params <- dbHatStar
  likeGrid_params$params <- params

  # Add level parameter corresponding to max likelihood.
  if(opt$levelParam) {
    # muHatStar  <- t(mu[, indexB, indexD, drop = FALSE])
    # mu is not a matrix so transpose t() doesnt apply.
    # It is a vector, so it should not matter much.
    # Also, without the drop = FALSE, it will automatically collapse to a vector.


    # # Troubleshooting here:
    # test <- array(seq(18), dim = c(2,3,3))
    # test
    # # Take a slice.
    # test_vec <- test[ ,2,2]
    # test_vec
    #
    # par_test <- matrix(c(100,200), nrow = 1, ncol = 2)
    # par_test
    #
    # c(par_test, test_vec)
    # cbind(par_test, test_vec)
    #
    # matrix(c(par_test, test_vec), nrow = 1, ncol = length(par_test) + length(test_vec))



    # print('mu[, indexB, indexD] = ')
    # print(mu[, indexB, indexD])
    # print('params = ')
    # print(params)


    muHatStar  <- mu[, indexB, indexD]
    likeGrid_params$muHatStar <- muHatStar


    # print('cbind(params, muHatStar) = ')
    # print(cbind(params, muHatStar))
    #
    # params <- cbind(params, muHatStar)

    # Don't ask:
    params <- matrix(c(params, muHatStar),
                     nrow = 1, ncol = length(params) + length(muHatStar))
    # Vectors and matrices and arrays. Oh my!


    likeGrid_params$params <- params
    # print('params = ')
    # print(params)

  }


  #--------------------------------------------------------------------------------
  # PLOT THE LIKELIHOOD
  #--------------------------------------------------------------------------------


  if (opt$plotLike) {

    plot.FCVARlikeGrid(likeGrid_params, k, r, opt, main = 'default')

  }


  return(likeGrid_params)
}


#' Plot the Likelihood Function for the FCVAR Model
#'
#' \code{plot.FCVARlikeGrid} plots the likelihood function from \code{FCVARlikeGrid}.
#' \code{FCVARlikeGrid} performs a grid-search optimization
#' by calculating the likelihood function
#' on a grid of candidate parameter values.
#' This function evaluates the likelihood over a grid of values
#' 	for \code{c(d,b)} (or \code{phi}).
#' 	It can be used when parameter estimates are sensitive to
#' 	starting values to give an approximation of the global max which can
#' 	then be used as the starting value in the numerical optimization in
#' 	\code{FCVARestn}.
#'
#' @param likeGrid_params A list output from \code{FCVARlikeGrid}.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @param file A string path and file name in which to plot a figure.
#' Renders to the plot window if running in RStudio if NULL.
#' @param file_ext A string file extension to indicate the graphics format.
#' Either png or pdf are currently supported.
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
#' opt$progress <- 2 # Show progress report on each value of b.
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, newOpt)
#' plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt, main = 'default')
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{plot.FCVARlikeGrid} plots the likelihood function from \code{FCVARlikeGrid}.
#' @export
#'
plot.FCVARlikeGrid <- function(likeGrid_params, k, r, opt,
                                file = NULL, file_ext = NULL,
                                main = NULL) {

  # Extract parameters.
  Grid2d <- likeGrid_params$Grid2d
  like <- likeGrid_params$like
  dGrid <- likeGrid_params$dGrid
  bGrid <- likeGrid_params$bGrid



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
      main <- c('Log-likelihood Function ',
                sprintf('Rank: %d, Lags: %d', r, k))
    }
  }


  # Plot likelihood depending on dimension of search.
  if(Grid2d) {
    # 2-dimensional plot.

    # Color palette (100 colors)
    # col.pal <- colorRampPalette(c("blue", "red"))
    # col.pal <- colorRampPalette(c("yellow", "red"))
    # colors <- col.pal(100)
    colors <- grDevices::rainbow(100)
    # colors <- heat.colors(100)
    # height of facets
    # like.facet.center <- (like[-1, -1] + like[-1, -ncol(like)] + like[-nrow(like), -1] + like[-nrow(like), -ncol(like)])/4
    like2D <- t(like)
    like.facet.center <- (like2D[-1, -1] + like2D[-1, -ncol(like2D)] + like2D[-nrow(like2D), -1] + like2D[-nrow(like2D), -ncol(like2D)])/4
    # like.facet.center <- like2D
    # Range of the facet center on a 100-scale (number of colors)
    like.facet.range <- cut(like.facet.center, 100)

    # Reduce the size of margins.
    # graphics::par()$mar
    # 5.1 4.1 4.1 2.1
    # bottom, left, top and right margins respectively.
    graphics::par(mar = c(1.1, 1.1, 2.1, 1.1))
    # Reset after.
    # graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

    graphics::persp(dGrid, bGrid,
          # like,
          like2D,
          xlab = 'd',
          ylab = 'b',
          # zlab = 'Log-likelihood',
          zlab = '',
          main = main,
          # phi = 45, theta = 45,
          phi = 45, theta = -55,
          # r = sqrt(3), # The distance of the eyepoint from the centre of the plotting box.
          r = 1.5,
          # d = 1, # The strength of the perspective transformation.
          d = 4.0,
          xaxs = "i",
          col = colors[like.facet.range]
    )

    # Reset the size of margins.
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

  } else {
    # 1-dimensional plot.
    if ((opt$R_psi[1] == 1) & (opt$R_psi[2] == -1) & (opt$r_psi[1] == 0)) {
      like_x_label <- 'd = b'
    } else {
      like_x_label <- 'phi'
    }

    graphics::plot(bGrid, like,
         main = main,
         ylab = 'log-likelihood',
         xlab = like_x_label,
         type = 'l',
         col = 'blue',
         lwd = 3)

  }

  # Close graphics device, if any.
  if (!is.null(file)) {
    grDevices::dev.off()
  }

}


#' Likelihood Function for the Unconstrained FCVAR Model
#'
#' \code{FCVARlikeMu} calculates the likelihood for the unconstrained FCVAR model
#' for a given set of parameter values. It is used by the \code{FCVARlikeGrid}
#' function to numerically optimize over the level parameter for given values of
#' the fractional parameters.
#'
#' @param y A matrix of variables to be included in the system.
#' @param db The orders of fractional integration.
#' @param mu The level parameter.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A number \code{like}, the log-likelihood evaluated at the
#' specified parameter values.
#' @examples
#' opt <- FCVARoptions()
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' like <- FCVARlikeMu(mu = colMeans(x), y = x, db = c(1, 1), k = 2, r = 1, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' The \code{FCVARlikeGrid} calls this function to perform a grid search over the
#' parameter values.
#' @export
#'
FCVARlikeMu <- function(mu, y, db, k, r, opt) {


  t <- nrow(y)
  x <- y - matrix(1, nrow = t, ncol = 1) %*% mu

  # print(summary(x))
  # print('mu = ')
  # print(mu)
  # print('db = ')
  # print(db)
  # print('k = ')
  # print(k)
  # print('r = ')
  # print(r)

  # Obtain concentrated parameter estimates.
  estimates <- GetParams(x, k, r, db, opt)

  # Calculate value of likelihood function.
  T <- t - opt$N
  p <- ncol(x)
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates$OmegaHat))



  return(like)
}


#' Likelihood Function for the Constrained FCVAR Model
#'
#' \code{FCVARlike} calculates the likelihood for the constrained FCVAR model
#' for a given set of parameter values. This function adjusts the variables with the level parameter,
#'if present, and returns the log-likelihood given \code{d} and \code{b}.
#'
#' @param params A vector of parameters \code{d} and \code{b}
#' (and \code{mu} if option selected).
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A number \code{like}, the log-likelihood evaluated at the
#' specified parameter values.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' FCVARlike(c(results$coeffs$db, results$coeffs$muHat), x, k = 2, r = 1, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' @export
#'
FCVARlike <- function(params, x, k, r, opt) {

  # global estimatesTEMP

  T <- nrow(x)
  p <- ncol(x)

  # print('params = ')
  # print(params)


  # If there is one linear restriction imposed, optimization is over phi,
  #  so translate from phi to (d,b), otherwise, take parameters as
  #  given.
  if(!is.null((opt$R_psi)) && nrow(opt$R_psi) == 1) {
    H <- pracma::null(opt$R_psi)
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

    # print('class(params) = ')
    # print(class(params))
    # print('class(mu) = ')
    # print(class(mu))
    #
    # print('mu = ')
    # print(mu)
    # print('diag(mu) = ')
    # print(diag(mu))
    #
    #
    # print('T = ')
    # print(T)
    y <- x - matrix(1, nrow = T, ncol = p) %*% diag(mu)
  } else {
    y <- x
  }


  # If d=b is not imposed via opt$restrictDB, but d>=b is required in
  #	opt$constrained, then impose the inequality here.
  if(opt$constrained) {
    b <- min(b, d)
  }


  # Obtain concentrated parameter estimates.
  estimates <- GetParams(y, k, r, c(d, b), opt)



  # Could have a function GetEstimates() that does all of the above.
  # estimates <- GetEstimates(params, x, k, r, opt)


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
  # T <- nrow(y) - opt$N
  # p <- ncol(y)
  # Base this on raw inputs:
  T <- nrow(x) - opt$N
  p <- ncol(x)
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates$OmegaHat))

  # Assign the value for the global variable estimatesTEMP.
  # Might adjust this later but it will work for now.
  # print('Search list in FCVARlike():')
  # print(search())
  # assign("estimatesTEMP", estimatesTEMP, envir = .GlobalEnv)

  return(like)
}


#' Obtain concentrated parameter estimates for the Constrained FCVAR Model
#'
#' \code{GetEstimates} calculates the concentrated parameter estimates for
#' the constrained FCVAR model for a given set of parameter values.
#' It is used after estimating the parameters that are numerically optimized
#' to obtain the corresponding concentrated parameters.
#' Like \code{FCVARlike} , this function adjusts the variables with the level parameter,
#' if present, and returns the log-likelihood given \code{d} and \code{b}.
#'
#' @param params A vector of parameters \code{d} and \code{b}
#' (and \code{mu} if option selected).
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A list object \code{estimates} containing the estimates,
#' including the following parameters:
#' \describe{
#'   \item{\code{db}}{The orders of fractional integration (taken directly from the input).}
#'   \item{\code{alphaHat}}{A \eqn{p x r} matrix of adjustment parameters.}
#'   \item{\code{betaHat}}{A \eqn{p x r} matrix of cointegrating vectors.
#'       The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.}
#'   \item{\code{rhoHat}}{A \eqn{p x 1} vector of restricted constatnts.}
#'   \item{\code{piHat}}{A \eqn{p x p} matrix \eqn{\Pi = \alpha \beta'} of long-run levels.}
#'   \item{\code{OmegaHat}}{A \eqn{p x p} covariance matrix of the error terms.}
#'   \item{\code{GammaHat}}{A ( \eqn{p x kp} matrix \code{cbind(GammaHat1,...,GammaHatk)})
#'   of autoregressive coefficients. }
#'   \item{\code{muHat}}{A vector of the optimal \code{mu} if level parameter is selected. }
#' }
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' GetEstimates(c(results$coeffs$db, results$coeffs$muHat), x, k = 2, r = 1, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARlike} performs the same calculations to obtain the value
#' of the likelihood function.
#' @export
#'
GetEstimates <- function(params, x, k, r, opt) {

  # global estimatesTEMP

  T <- nrow(x)
  p <- ncol(x)

  # print('params = ')
  # print(params)


  # If there is one linear restriction imposed, optimization is over phi,
  #  so translate from phi to (d,b), otherwise, take parameters as
  #  given.
  if(!is.null((opt$R_psi)) && nrow(opt$R_psi) == 1) {
    H <- pracma::null(opt$R_psi)
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

    # print('class(params) = ')
    # print(class(params))
    # print('class(mu) = ')
    # print(class(mu))
    #
    # print('mu = ')
    # print(mu)
    # print('diag(mu) = ')
    # print(diag(mu))
    #
    #
    # print('T = ')
    # print(T)
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
  # estimatesTEMP <- estimates
  # If level parameter is present, they are the last p parameters in the
  #   params vector
  if (opt$levelParam) {
    estimates$muHat <- mu
  } else {
    estimates$muHat <- NULL
  }

  return(estimates)

}


#' Likelihood Function for the FCVAR Model
#'
#' \code{FCVARlikeFull} calculates the likelihood for the constrained FCVAR model
#' for a given set of parameter values.
#' This function takes the full set of coefficients from a call
#' to \code{FCVARestn} and returns the log-likelihood given \code{d} and \code{b}.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param coeffs A list of coefficients of the FCVAR model.
#' It is an element of the list \code{results} returned by \code{FCVARestn},
#' without the parameters \code{betaHat} and \code{rhoHat}.
#' @param betaHat A \eqn{p x r} matrix of cointegrating vectors.
#' The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.
#' @param rhoHat A \eqn{p x 1} vector of restricted constatnts.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A number \code{like}, the log-likelihood evaluated at the
#' specified parameter values.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' FCVARlikeFull(x, k = 2, r = 1, coeffs = results$coeffs,
#'               beta = results$coeffs$betaHat, rho = results$coeffs$rhoHat, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} for the estimation of coefficients in \code{coeffs}.
#' @export
#'
FCVARlikeFull <- function(x, k, r, coeffs, betaHat, rhoHat, opt) {


  # Add betaHat and rhoHat to the coefficients to get residuals because they
  #   are not used in the Hessian calculation and are missing from the
  #   structure coeffs
  coeffs$betaHat <- betaHat
  coeffs$rhoHat <- rhoHat

  # Obtain residuals.
  epsilon <- GetResiduals(x, k, r, coeffs, opt)


  # Calculate value of likelihood function.
  T <- nrow(x) - opt$N
  p <- ncol(x)
  OmegaHat <- t(epsilon) %*% epsilon/T
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(OmegaHat))

  return(like)
}


#' Transform Data for Regression
#'
#' \code{TransformData} transforms the raw data by fractional differencing.
#' The output is in the format required for regression and
#' reduced rank regression.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param db The orders of fractional integration.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A list object \code{Z_array} containing the transformed data,
#' including the following parameters:
#' \describe{
#'   \item{\code{Z0}}{A matrix of data calculated by fractionally differencing
#'   \code{x} at differencing order \code{d}.}
#'   \item{\code{Z1}}{The matrix \code{x} (augmented with a vector of ones
#'   if the model includes a restricted constant term), which is then lagged and stacked
#'   and fractionally differenced (with order \eqn{d-b}).}
#'   \item{\code{Z2}}{The matrix \code{x} lagged and stacked
#'   and fractionally differenced (with order \code{d}).}
#'   \item{\code{Z3}}{A column of ones if model includes an
#'   unrestricted constant term, otherwise \code{NULL}.}
#' }
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' Z_array <- TransformData(x, k = 2, db = results$coeffs$db, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} calls \code{GetParams}, which calls \code{TransformData}
#' to estimate the FCVAR model.
#' \code{TransformData} in turn calls \code{FracDiff} and \code{Lbk}
#' to perform the transformation.
#' @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
#' "Likelihood inference for a fractionally cointegrated
#' vector autoregressive model," Econometrica 80, 2667-2732.
#' @references Johansen, S. (1995). "Likelihood-Based Inference
#' in Cointegrated Vector Autoregressive Models,"
#' New York: Oxford University Press.
#' @export
#'
TransformData <- function(x, k, db, opt) {


  # print(summary(x))

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


#' Calculate Residuals for the FCVAR Model
#'
#' \code{GetResiduals} calculates residuals for the FCVAR model
#' from given parameter values.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param coeffs A list of coefficients of the FCVAR model.
#' It is an element of the list \code{results} returned by \code{FCVARestn}.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A matrix \code{epsilon} of residuals from FCVAR model estimation
#' calculated with the parameter estimates specified in \code{coeffs}.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' epsilon <- GetResiduals(x, k = 2, r = 1, coeffs = results$coeffs, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} to estimate the FCVAR model.
#' @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
#' "Likelihood inference for a fractionally cointegrated
#' vector autoregressive model," Econometrica 80, 2667-2732.
#' @export
#'
GetResiduals <- function(x, k, r, coeffs, opt) {


  # If level parameter is included, the data must be shifted before
  #   calculating the residuals:
  if (opt$levelParam) {
    T <- nrow(x)
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
      Z1 %*% rbind(coeffs$betaHat, coeffs$rhoHat) %*% t(coeffs$alphaHat)
  }


  if (k > 0) {
    epsilon <- epsilon - Z2 %*% t(coeffs$GammaHat)
  }


  if (opt$unrConstant) {
    epsilon <- epsilon - Z3 %*% t(coeffs$xiHat)
  }



  return(epsilon)
}


#' Calculate Lag Polynomial in the Fractional Lag Operator
#'
#' \code{Lbk} calculates a lag polynomial in the fractional lag operator.
#'
#' @param x A matrix of variables to be included in the system.
#' @param b The order of fractional differencing.
#' @param k The number of lags in the system.
#' @return A matrix \code{Lbkx} of the form \eqn{[ Lb^1 x, Lb^2 x, ..., Lb^k x]}
#' where \eqn{Lb = 1 - (1-L)^b}.
#' The output matrix has the same number of rows as \code{x}
#' but \code{k} times as many columns.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' Lbkx <- Lbk(x, b = results$coeffs$db[2], k = 2)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} calls \code{GetParams}, which calls \code{TransformData}
#' to estimate the FCVAR model.
#' \code{TransformData} in turn calls \code{FracDiff} and \code{Lbk}
#' to perform the transformation.
#' @export
#'
Lbk <- function(x, b, k) {

  # print('summary(x) = ')
  # print(summary(x))
  # print('b = ')
  # print(b)
  # print('k = ')
  # print(k)

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
                    ( Lbkx[ , (p*(i-2)+1) : ncol(Lbkx)] -
                        FracDiff(Lbkx[ , (p*(i-2)+1) : ncol(Lbkx)], b) ))

    }

  }



  return(Lbkx)
}



#' Fast Fractional Differencing
#'
#' \code{FracDiff} is a fractional differencing procedure based on the
#' 	fast fractional difference algorithm of Jensen & Nielsen (2014).
#'
#' @param x A matrix of variables to be included in the system.
#' @param d The order of fractional differencing.
#' @return A vector or matrix \code{dx} equal to \eqn{(1-L)^d x}
#' of the same dimensions as x.
#' @examples
#' set.seed(42)
#' WN <- matrix(stats::rnorm(200), nrow = 100, ncol = 2)
#' \dontrun{MVWNtest_stats <- MVWNtest(x = WN, maxlag = 10, printResults = 1)}
#' x <- FracDiff(x = WN, d = - 0.5)
#' \dontrun{MVWNtest_stats <- MVWNtest(x = x, maxlag = 10, printResults = 1)}
#' WN_x_d <- FracDiff(x, d = 0.5)
#' \dontrun{MVWNtest_stats <- MVWNtest(x = WN_x_d, maxlag = 10, printResults = 1)}
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} calls \code{GetParams}, which calls \code{TransformData}
#' to estimate the FCVAR model.
#' \code{TransformData} in turn calls \code{FracDiff} and \code{Lbk}
#' to perform the transformation.
#' @references Jensen, A. N. and M. \enc{Ø}{O}. Nielsen (2014).
#' "A fast fractional difference algorithm,"
#' Journal of Time Series Analysis 35, 428-436.
#' @note This function differs from the \code{diffseries} function
#' in the \code{fracdiff} package, in that the \code{diffseries}
#' function demeans the series first.
#' In particular, the function calls \code{FCVAR::FracDiff(x - mean(x), d = 0.5)}
#' and \code{fracdiff::diffseries(x, d = 0.5)} return the same values.
#' (NOTE TO US: We should probably reconsider our naming choice for this function,
#' since the \code{fracdiff} function in the \code{fracdiff} package
#' actually estimates the ARFIMA model. Not an ideal choice of name
#' on their side but they came first.)
#' @export
#'
FracDiff <- function(x, d) {


  # print(summary(x))

  if(is.null(x)) {
    dx <- NULL
  } else {


    # T <- nrow(x)
    p <- ncol(x)



    iT <- nrow(x)
    np2 <- stats::nextn(2*iT - 1, 2)

    k <- 1:(iT-1)
    b <- c(1, cumprod((k - d - 1)/k))


    # print('np2 = ')
    # print(np2)
    # print('iT = ')
    # print(iT)
    #
    # print('c(...) = ')
    # print(summary(c(x, rep(0, np2 - iT))))

    dx <- matrix(0, nrow = iT, ncol = p)

    for (i in 1:p) {


      dxi <- stats::fft(stats::fft(c(b, rep(0, np2 - iT))) *
                          stats::fft(c(x[, i], rep(0, np2 - iT))), inverse = T) / np2

      dx[, i] <- Re(dxi[1:iT])

    }



  }

  return(dx)
}


#' Count the Number of Free Parameters
#'
#' \code{GetFreeParams} counts the number of free parameters based on
#' 	the number of coefficients to estimate minus the total number of
#' 	restrictions. When both \code{alpha} and \code{beta} are restricted,
#' 	the rank condition is used to count the free parameters in those two variables.
#'
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param p The number of variables in the system.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @param rankJ The rank of a conditioning matrix, as described in
#' Boswijk & Doornik (2004, p.447), which is only used if there are
#' restrictions imposed on \code{alpha} or \code{beta}, otherwise \code{NULL}.
#' @return The number of free parameters \code{fp}.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = NULL)
#'
#' opt <- FCVARoptions()
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = 4)
#'
#' opt <- FCVARoptions()
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
#' newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = 4)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn}, \code{HypoTest} and \code{LagSelect} to estimate the FCVAR model
#' and use this in the calculation of the degrees of freedom
#' for a variety of statistics.
#' @references Boswijk, H. P. and J. A. Doornik (2004).
#' "Identifying, estimating and testing restricted cointegrated systems:
#' An overview," Statistica Neerlandica 58, 440-465.
#' @export
#'
GetFreeParams <- function(k, r, p, opt, rankJ) {

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


  if(is.null(opt$R_Beta) & is.null(opt$R_Alpha)) {
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


#' Calculate the Hessian Matrix
#'
#' \code{FCVARhess} calculates the Hessian matrix of the
#' 	log-likelihood by taking numerical derivatives.
#' 	It is used to calculate the standard errors of parameter estimates.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param coeffs A list of coefficients of the FCVAR model.
#' It is an element of the list \code{results} returned by \code{FCVARestn}.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return The \code{hessian} matrix  of second derivatives of the FCVAR
#' log-likelihood function, calculated with the parameter estimates
#' specified in \code{coeffs}.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' hessian <- FCVARhess(x, k = 2, r = 1, coeffs = results$coeffs, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} to estimate the FCVAR model and calculate
#' standard errors of the estimates.
#' @export
#'
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
      like1 <- FCVARlikeFull(x, k, r, coeffsAdj, beta, rho, opt )

      # negative shift in first parameter, positive shift in second
      # parameter. If same parameter, no shift.
      phi2 <- phi0 - deltaPhi*( (1:nPhi) == i ) + deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi2, k, r, p, opt )
      # calculate likelihood
      like2 <- FCVARlikeFull(x, k, r, coeffsAdj, beta, rho, opt)

      # positive shift in first parameter, negative shift in second
      # parameter. If same parameter, no shift.
      phi3 <- phi0 + deltaPhi*( (1:nPhi) == i ) - deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi3, k, r, p, opt )
      # calculate likelihood
      like3 <- FCVARlikeFull(x, k, r, coeffsAdj, beta, rho, opt)

      # negative shift in both parameters.
      phi4 <- phi0 - deltaPhi*( (1:nPhi) == i ) - deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi4, k, r, p, opt )
      # calculate likelihood
      like4 <- FCVARlikeFull(x, k, r, coeffsAdj, beta, rho, opt)

      # The numerical approximation to the second derivative.
      hessian[i,j] <- ( like1 - like2 - like3 + like4 )/4/deltaPhi[i]/deltaPhi[j]

      # Hessian is symmetric.
      hessian[j,i] <- hessian[i,j]


    }

  }

  return(hessian)
}


#' Collect Parameters into a Vector
#'
#' \code{SEmat2vecU} transforms the model parameters in matrix
#' 	form into a vector.
#'
#' @param coeffs A list of coefficients of the FCVAR model.
#' It is an element of the list \code{results} returned by \code{FCVARestn}.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param p The number of variables in the system.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return A vector \code{param} of parameters in the FCVAR model.
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' params <- SEmat2vecU(coeffs = results$coeffs, k = 2, r = 1, p = 3, opt)
#' coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
#'
#' params <- matrix(seq(25))
#' coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
#' params <- SEmat2vecU(coeffs = coeffs, k = 2, r = 1, p = 3, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} to estimate the FCVAR model and calculate
#' standard errors of the estimates.
#' \code{SEmat2vecU} is called by \code{FCVARhess} to sort the parameters
#' into a vector to calculate the Hessian matrix.
#' \code{SEvec2matU} is a near inverse of \code{SEmat2vecU},
#' in the sense that \code{SEvec2matU} obtains only a
#' subset of the parameters in \code{results$coeffs}.
#' @export
#'
SEmat2vecU <- function(coeffs, k, r, p , opt) {

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



#' Extract Parameters from a Vector
#'
#' \code{SEvec2matU} transforms the vectorized model parameters
#' 	into matrices.
#'
#' @param param A vector of parameters in the FCVAR model.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param p The number of variables in the system.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return \code{coeffs}, a list of coefficients of the FCVAR model.
#' It has the same form as an element of the list \code{results}
#' returned by \code{FCVARestn}.
#'
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' params <- SEmat2vecU(coeffs = results$coeffs, k = 2, r = 1, p = 3, opt)
#' coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
#'
#' params <- matrix(seq(25))
#' coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
#' params <- SEmat2vecU(coeffs = coeffs, k = 2, r = 1, p = 3, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} to estimate the FCVAR model and calculate
#' standard errors of the estimates.
#' \code{SEmat2vecU} is called by \code{FCVARhess} to convert the parameters
#' from a vector into the coefficients after calculating the Hessian matrix.
#' \code{SEmat2vecU} is a near inverse of \code{SEvec2matU},
#' in the sense that \code{SEvec2matU} obtains only a
#' subset of the parameters in \code{results$coeffs}.
#' @export
#'
SEvec2matU <- function(param, k, r, p, opt ) {

  # Create list for output.
  coeffs <- list(
    db = NULL,
    muHat = NULL,
    xiHat = NULL,
    alphaHat = NULL,
    GammaHat = NULL
  )

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




#' Calculate Restricted Estimates for the FCVAR Model
#'
#' \code{GetRestrictedParams} calculates restricted estimates of
#' cointegration parameters and error variance in the FCVAR model.
#' To calculate the optimum, it uses the switching algorithm
#' of Boswijk and Doornik (2004, page 455) to optimize over free parameters
#' \eqn{\psi} and \eqn{\phi} directly, combined with the line search proposed by
#' Doornik (2016, working paper). We translate between  \eqn{(\psi, \phi)} and
#' \eqn{(\alpha, \beta)} using the relation of \eqn{R_\alpha vec(\alpha) = 0} and
#' \eqn{A\psi = vec(\alpha')}, and \eqn{R_\beta vec(\beta) = r_\beta} and
#' \eqn{H\phi + h = vec(\beta)}. Note the transposes.
#'
#' @param beta0 The unrestricted estimate of \code{beta},
#' a \eqn{p x r} matrix of cointegrating vectors
#' returned from \code{FCVARestn} or \code{GetParams}.
#' @param S00 A matrix of product moments,
#' calculated from the output of \code{TransformData} in \code{GetParams}.
#' @param S01 A matrix of product moments,
#' calculated from the output of \code{TransformData} in \code{GetParams}.
#' @param S11 A matrix of product moments,
#' calculated from the output of \code{TransformData} in \code{GetParams}.
#' @param T The number of observations in the sample.
#' @param p The number of variables in the system.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#'
#'#' @return A list object \code{switched_mats} containing the restricted estimates,
#' including the following parameters:
#' \describe{
#'   \item{\code{betaStar}}{A \eqn{p x r} matrix of cointegrating vectors.
#'       The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.}
#'   \item{\code{alphaHat}}{A \eqn{p x r} matrix of adjustment parameters.}
#'   \item{\code{OmegaHat}}{A \eqn{p x p} covariance matrix of the error terms.}
#' }
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' opt$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
#' opt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' betaStar <- matrix(c(-0.95335616, -0.07345676, -0.29277318), nrow = 3)
#' S00 <- matrix(c(0.0302086527,  0.001308664,  0.0008200923,
#'                 0.0013086640,  0.821417610, -0.1104617893,
#'                 0.0008200923, -0.110461789,  0.0272861128), nrow = 3)
#' S01 <- matrix(c(-0.0047314320, -0.04488533,  0.006336798,
#'                 0.0026708007,  0.17463884, -0.069006455,
#'                 -0.0003414163, -0.07110324,  0.022830494), nrow = 3, byrow = TRUE)
#' S11 <- matrix(c( 0.061355941,  -0.4109969,  -0.007468716,
#'                  -0.410996895,  70.6110313, -15.865097810,
#'                  -0.007468716, -15.8650978,   3.992435799), nrow = 3)
#' switched_mats <- GetRestrictedParams(betaStar, S00, S01, S11, T = 316, p = 3, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{FCVARestn} calls \code{GetParams} to estimate the FCVAR model,
#' which in turn calls \code{GetRestrictedParams} if there are restrictions
#' imposed on \code{alpha} or \code{beta}.
#' @references Boswijk, H. P. and J. A. Doornik (2004).
#' "Identifying, estimating and testing restricted cointegrated systems:
#' An overview," Statistica Neerlandica 58, 440-465.
#' @references Doornik, J. A. (2018).
#' "Accelerated estimation of switching algorithms: the cointegrated
#' VAR model and other applications,"
#' Forthcoming in Scandinavian Journal of Statistics.
#' @export
#'
GetRestrictedParams <- function(beta0, S00, S01, S11, T, p, opt) {


  r  <- ncol(beta0)
  p1 <- p + opt$rConstant

  # Restrictions on beta.
  if(is.null(opt$R_Beta)) {
    H <- diag(p1*r)
    h <- matrix(0, nrow = p1*r, ncol = 1)
  } else {
    H <- pracma::null(opt$R_Beta)
    h <- t(opt$R_Beta) %*%
      solve(opt$R_Beta %*% t(opt$R_Beta)) %*% opt$r_Beta
  }


  # Restrictions on alpha.
  #   We use the commutation matrix K_pr to transform vec(A) into vec(A'),
  #   see Magnus & Neudecker (1988, p. 47, eqn (1)).
  Ip  <- diag(p)
  Kpr <- matrix(kronecker(Ip, diag(r)), nrow = p*r, ncol = p*r)
  if(is.null(opt$R_Alpha)) {
    A <- Kpr %*% diag(p*r)
  } else {
    A <- pracma::null(opt$R_Alpha %*% solve(Kpr))
  }


  # Least squares estimator of Pi, used in calculations below.
  PiLS <- t((S01 %*% solve(S11)))
  vecPiLS <- matrix(PiLS, nrow = length(PiLS), ncol = 1)





  # Starting values for switching algorithm.
  betaStar <- beta0
  alphaHat <-  S01 %*% beta0 %*% solve(t(beta0) %*% S11 %*% beta0)
  OmegaHat <- S00 - S01 %*% betaStar %*% t(alphaHat) -
    alphaHat %*% t(betaStar) %*% t(S01) +
    alphaHat %*% t(betaStar) %*% S11 %*% betaStar %*% t(alphaHat)

  # Algorithm specifications.
  iters   <- opt$UncFminOptions$MaxFunEvals
  Tol     <- opt$UncFminOptions$TolFun
  conv    <- 0
  i       <- 0

  # Tolerance for entering the line search epsilon_s in Doornik's paper.
  TolSearch <- 0.01

  # Line search parameters.
  lambda <- c(1, 1.2, 2, 4, 8)
  nS <- length(lambda)
  likeSearch <- matrix(NA, nrow = nS, ncol = 1)
  OmegaSearch <- array(0, dim = c(p, p,nS))

  # print('vecPiLS = ')
  # print(vecPiLS)
  # print('alphaHat = ')
  # print(alphaHat)
  # print('kronecker(alphaHat, diag(p1)) = ')
  # print(kronecker(alphaHat, diag(p1)))
  # print('h = ')
  # print(h)
  # print('kronecker(alphaHat, diag(p1)) %*% h = ')
  # print(kronecker(alphaHat, diag(p1)) %*% h)
  #



  # Get candidate values for entering the switching algorithm.
  vecPhi1 <- solve(t(H) %*%
                     kronecker(t(alphaHat) %*% solve(OmegaHat) %*% alphaHat, S11) %*%
                     H) %*%
    t(H) %*% (kronecker(t(alphaHat) %*% solve(OmegaHat), S11)) %*%
    (vecPiLS - kronecker(alphaHat, diag(p1)) %*% h)

  # Translate vecPhi to betaStar.
  vecB <- H %*% vecPhi1 + h
  betaStar <- matrix(vecB, nrow = p1, ncol = r)

  # Candidate value of vecPsi.
  vecPsi1 <- solve(t(A) %*%
                     kronecker(solve(OmegaHat), t(betaStar) %*% S11 %*% betaStar) %*%
                     A) %*%
    t(A) %*% (kronecker(solve(OmegaHat), t(betaStar) %*% S11)) %*% vecPiLS

  # Translate vecPsi to alphaHat.
  vecA <- A %*% vecPsi1 # This is vec(alpha')
  alphaHat <- matrix(solve(Kpr) %*% vecA, nrow = p, ncol = r)

  # Candidate values of piHat and OmegaHat.
  piHat1 <- alphaHat %*% t(betaStar)
  OmegaHat <- S00 - S01 %*% betaStar %*% t(alphaHat) -
    alphaHat %*% t(betaStar) %*% t(S01) +
    alphaHat %*% t(betaStar) %*% S11 %*% betaStar %*% t(alphaHat)





  # Calculate the likelihood.
  like1 <- - log(det(OmegaHat))

  while(i<=iters && !conv) {

    # Update values for convergence criteria.
    piHat0  <- piHat1
    like0   <- like1

    if(i == 1) {
      # Initialize candidate values.
      vecPhi0_c <- vecPhi1
      vecPsi0_c <- vecPsi1
    }


    #  ---- alpha update step ---- %
    # Update vecPsi.
    vecPsi1 <- solve(t(A) %*% kronecker(solve(OmegaHat), t(betaStar) %*%
                                     S11 %*% betaStar) %*% A) %*%
      t(A) %*% (kronecker(solve(OmegaHat), t(betaStar) %*% S11)) %*% vecPiLS

    # Translate vecPsi to alphaHat.
    vecA <- A %*% vecPsi1 # This is vec(alpha')
    alphaHat <- matrix(solve(Kpr) %*% vecA, nrow = p, ncol = r)

    #  ---- omega update step ---- %
    # Update OmegaHat.
    OmegaHat <- S00 - S01 %*% betaStar %*% t(alphaHat) -
      alphaHat %*% t(betaStar) %*% t(S01) +
      alphaHat %*% t(betaStar) %*% S11 %*% betaStar %*% t(alphaHat)

    #  ---- beta update step ---- %

    # Update vecPhi.
    vecPhi1 <- solve(t(H) %*%
                       kronecker(t(alphaHat) %*% solve(OmegaHat) %*% alphaHat, S11) %*%
                       H) %*%
      t(H) %*% (kronecker(t(alphaHat) %*% solve(OmegaHat), S11)) %*%
      (vecPiLS - kronecker(alphaHat, diag(p1)) %*% h)

    # Translate vecPhi to betaStar.
    vecB <- H %*% vecPhi1 + h
    betaStar <- matrix(vecB, nrow = p1, ncol = r)

    #  ---- pi and likelihood update  ---- %
    # Update estimate of piHat.
    piHat1 <- alphaHat %*% t(betaStar)

    # Update OmegaHat with new alpha and new beta.
    OmegaHat <- S00 - S01 %*% betaStar %*% t(alphaHat) -
      alphaHat %*% t(betaStar) %*% t(S01) +
      alphaHat %*% t(betaStar) %*% S11 %*% betaStar %*% t(alphaHat)




    # Calculate the likelihood.
    like1 <- - log(det(OmegaHat))

    if(i > 0) {


      # Calculate relative change in likelihood
      likeChange <- (like1 - like0) / (1 + abs(like0))



      # Check relative change and enter line search if below tolerance.
      if(likeChange < TolSearch && opt$LineSearch) {


        # Calculate changes in parameters.
        deltaPhi <- vecPhi1 - vecPhi0_c
        deltaPsi <- vecPsi1 - vecPsi0_c

        # Initialize parameter bins
        vecPhi2 <- matrix(NA, nrow = length(vecPhi1), ncol = nS)
        vecPsi2 <- matrix(NA, nrow = length(vecPsi1), ncol = nS)

        # Values already calculated for lambda = 1.
        vecPhi2[ ,1]       <- vecPhi1
        vecPsi2[ ,1]       <- vecPsi1
        likeSearch[1]      <- like1
        OmegaSearch[ , ,1] <- OmegaHat



        for (iL in 2:nS) {

          # New candidates for parameters based on line search.
          vecPhi2[ , iL] <- vecPhi0_c + lambda[iL]*deltaPhi
          vecPsi2[ , iL] <- vecPsi0_c + lambda[iL]*deltaPsi

          # Translate to alpha and beta
          vecA <- A %*% vecPsi2[ , iL]
          alphaHat <- matrix(solve(Kpr) %*% vecA, nrow = p, ncol = r)
          vecB <- H %*% vecPhi2[ , iL] + h
          betaStar <- matrix(vecB, nrow = p1, ncol = r)

          # Calculate and store OmegaHat
          OmegaSearch[ , , iL] <- S00 -
            S01 %*% betaStar %*% t(alphaHat) -

            alphaHat %*% t(betaStar) %*% t(S01) +
            alphaHat %*% t(betaStar) %*% S11 %*% betaStar %*% t(alphaHat)

          # Calculate and store log-likelihood
          likeSearch[iL] <- - log(det(OmegaSearch[ , , iL]))

        }




        # Update max likelihood and OmegaHat based on line search.
        # [ iSearch ] <- find(likeSearch == max(max(likeSearch)))
        iSearch <- which(likeSearch == max(likeSearch), arr.ind = TRUE)

        # If there are identical likelihoods, choose smallest
        # increment
        iSearch  <- min(iSearch)
        like1    <- likeSearch[iSearch]
        OmegaHat <- OmegaSearch[ , , iSearch]
        # Save old candidate parameter vectors for next iteration.
        vecPhi0_c  <- vecPhi1
        vecPsi0_c  <- vecPsi1

        # Update new candidate parameter vectors.
        vecPhi1  <- vecPhi2[ ,iSearch]
        vecPsi1  <- vecPsi2[ ,iSearch]
        # Update coefficients.
        vecA     <- A %*% vecPsi1
        alphaHat <- matrix(solve(Kpr) %*% vecA, nrow = p, ncol = r)
        vecB     <- H %*% vecPhi1 + h
        betaStar <- matrix(vecB, nrow = p1, ncol = r)
        # Update estimate of piHat.
        piHat1 <- alphaHat %*% t(betaStar)


      }


      # Calculate relative change in likelihood
      likeChange <- (like1 - like0) / (1 + abs(like0))

      # Calculate relative change in coefficients.
      piChange   <- max(abs(piHat1 - piHat0) / (1 + abs(piHat0)) )

      # Check convergence.
      if(abs(likeChange) <= Tol && piChange <= sqrt(Tol)) {
        conv <- 1
      }


    }

    i <- i+1
  }




  # Return a list of parameters from the switching algorithm.
  switched_mats <- list(
    betaStar = betaStar,
    alphaHat = alphaHat,
    OmegaHat = OmegaHat
  )

  return(switched_mats)
}






