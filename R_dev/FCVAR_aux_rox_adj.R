#' Draw Samples from the FCVAR Model
#'
#' \code{FCVARsim} simulates the FCVAR model as specified by
#' input \code{model} and starting values specified by \code{data}.
#' Errors are drawn from a normal distribution.
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
#' It creates a wild bootstrap sample by augmenting each iteration
#' with a bootstrap error. The errors are sampled from the
#' residuals specified under the \code{model} input and have a
#' positive or negative sign with equal probability (the Rademacher distribution).
#' @param data A \eqn{T x p} matrix of starting values for the simulated realizations.
#' @param model A list of estimation results, just as if estimated from \code{FCVARest}.
#' The parameters in \code{model} can also be set or adjusted by assigning new values.
#' @param NumPeriods The number of time periods in the simulation.
#' @return A \code{NumPeriods} by \eqn{p} matrix \code{xBS} of simulated bootstrap values.
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

  # Generate disturbance term for use in the bootstrap.

  # Center residuals
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


    # Append x with zeros to simplify calculations
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


# Calculate Parameters of the Cointegrated VAR Model
#
# \code{GetParams} calculates the parameters of the cointegrated VAR model
#  from a multivariate sample.
#  This function uses FWL and reduced rank regression to obtain
# the estimates of \code{Alpha}, \code{Beta}, \code{Rho}, \code{Pi},
# \code{Gamma}, and \code{Omega}, using the notation in Johansen and Nielsen (2012).
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param x A matrix of variables to be included in the system.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param db The orders of fractional integration.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A list object \code{estimates} containing the estimates,
# including the following parameters:
# \describe{
#   \item{\code{db}}{The orders of fractional integration (taken directly from the input).}
#   \item{\code{alphaHat}}{A \eqn{p x r} matrix of adjustment parameters.}
#   \item{\code{betaHat}}{A \eqn{p x r} matrix of cointegrating vectors.
#       The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.}
#   \item{\code{rhoHat}}{A \eqn{p x 1} vector of restricted constatnts.}
#   \item{\code{piHat}}{A \eqn{p x p} matrix \eqn{\Pi = \alpha \beta'} of long-run levels.}
#   \item{\code{OmegaHat}}{A \eqn{p x p} covariance matrix of the error terms.}
#   \item{\code{GammaHat}}{A ( \eqn{p x kp} matrix \code{cbind(GammaHat1,...,GammaHatk)})
#   of autoregressive coefficients. }
# }
# @examples
# opt <- FCVARoptions()
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# estimates <- GetParams(x, k = 2, r = 1, db = c(1, 1), opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} calls \code{GetParams} to estimate the FCVAR model.
# @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
# "Likelihood inference for a fractionally cointegrated
# vector autoregressive model," Econometrica 80, 2667-2732.
# @export
#
GetParams <- function(x, k, r, db, opt) {


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

  if ( (k == 0) || is.infinite(kappa(t(Z2hat) %*% Z2hat)) ) {

    # No lags, so there are no effects of Z2.
    R0 <- Z0hat
    R1 <- Z1hat
    # Includes the case when Z2hat is near zero and
    # t(Z2hat) %*% Z2hat is ill-conditioned.

  } else {

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
    betaStar <- NA
    alphaHat <- NA
    PiHat <- NULL
    rhoHat <- NULL
    OmegaHat <- S00

  } else if (( r > 0 ) & ( r < p )) {


    # Calculate solution from eigenvectors
    eig_out <- eigen( solve(S11) %*% S10 %*% solve(S00) %*% S01 )
    D <- eig_out$values # Only the vector, not the diagonal matrix.
    V1 <- eig_out$vectors

    V2 <- t( cbind(t(V1), D)[order(D), ] )

    V <- V2[1:p1, , drop = FALSE]
    betaStar <- V[ 1:p1, seq(p1, p1-r+1, by = -1), drop = FALSE]



    # If either alpha or beta is restricted, then the likelihood is
    #   maximized subject to the constraints. This section of the code
    #   replaces the call to the switching algorithm in the previous
    #   version.
    if (!is.null(opt$R_Alpha) | !is.null(opt$R_Beta) ) {


      switched_mats <- GetRestrictedParams(betaStar, S00, S01, S11, T, p, opt)
      betaStar <- switched_mats$betaStar
      alphaHat <- switched_mats$alphaHat
      OmegaHat <- switched_mats$OmegaHat

      betaHat <- betaStar[1:p, 1:r, drop = FALSE]
      PiHat <- alphaHat %*% t(betaHat)



    } else {

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

    # Check conformability:
    betaStar <- rbind(betaHat, rhoHat)

  }

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
#' 	\code{plot.FCVAR_grid} plots the likelihood function from \code{FCVARlikeGrid}.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
#' generated from \code{FCVARoptions()}.
#' @return An S3 object of type \code{FCVAR_grid} containing the optimization results,
#' including the following parameters:
#' \describe{
#'   \item{\code{params}}{A vector \code{params} of \code{d} and \code{b}
#'   (and \code{mu} if level parameter is selected)
#'   corresponding to a maximum over the grid of \code{c(d,b)} or \code{phi}.}
#'   \item{\code{dbHatStar}}{A vector of \code{d} and \code{b}
#'   corresponding to a maximum over the grid of \code{c(d,b)} or \code{phi}.}
#'   \item{\code{muHatStar}}{A vector of the optimal \code{mu} if level parameter is selected. }
#'   \item{\code{Grid2d}}{An indicator for whether or not the optimization
#'   is conducted over a 2-dimensional parameter space,
#'   i.e. if there is no equality restriction on \code{d} and \code{b}.}
#'   \item{\code{dGrid}}{A vector of the grid points in the parameter \code{d},
#'     after any transformations for restrictions, if any.}
#'   \item{\code{bGrid}}{A vector of the grid points in the parameter \code{b},
#'     after any transformations for restrictions, if any.}
#'   \item{\code{dGrid_orig}}{A vector of the grid points in the parameter \code{d},
#'     in units of the fractional integration parameter.}
#'   \item{\code{bGrid_orig}}{A vector of the grid points in the parameter \code{b},
#'     in units of the fractional integration parameter.}
#'   \item{\code{like}}{The maximum value of the likelihood function over the chosen grid.}
#'   \item{\code{k}}{The number of lags in the system.}
#'   \item{\code{r}}{The cointegrating rank.}
#'   \item{\code{opt}}{An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
#'     generated from \code{FCVARoptions()}.}
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
#' # newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, opt)
#' \dontrun{plot.FCVAR_grid(likeGrid_params)}
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
#' # newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, opt)
#' \dontrun{plot.FCVAR_grid(likeGrid_params)}
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
#' # newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, opt)
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
#' # newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, opt)
#' }
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{plot.FCVAR_grid} plots the likelihood function from \code{FCVARlikeGrid}.
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


  #--------------------------------------------------------------------------------
  # INITIALIZATION
  #--------------------------------------------------------------------------------

  # Check if performing a 2-dimensional or 1-dimensional grid search. The
  #   step size will be smaller if searching in only one dimension.
  if(is.null(opt$R_psi)) {
    Grid2d <- 1
    dbStep <- opt$dbStep2D # This variable was added to options.
  } else {
    Grid2d <- 0
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
      # following block calculates the total number of iterations
      # by counting the number of instances in which the grid values
      # in b are less than or equal to those in d. This should be
      # moved to the beginning of the loop for faster computation
      # later.
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

  # Initialize progress reports.
  if(opt$progress != 0) {
    if(opt$progress == 1) {
      progbar <- utils::txtProgressBar(min = 0, max = totIters, style = 3)

    }

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

    # Print a progress report of the preferred style.
    if (opt$progress == 1) {
      utils::setTxtProgressBar(progbar, value = iterCount)
    } else if (opt$progress == 2) {
      message(sprintf('Now estimating for iteration %d of %d: b = %f.',
                      iterCount, totIters, b))
    }

    # If d>=b (constrained) then search only over values of d that are
    #   greater than or equal to b. Also, perform this operation only if
    #   searching over both d and b.
    if(opt$constrained & Grid2d) {
      # Start at the index that corresponds to the first value in the
      # grid for d that is >= b
      dStart <- which(dGrid >= b)[1]
    }


    for (iD in dStart:nD) {

      iterCount <- iterCount + 1

      # db is defined in terms of d,b for use by FCVARlikeMU
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

        min_out <- stats::optim(StartVal,
                         function(params) {-FCVARlikeMu(params, x, db, k, r, opt)})

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


    }


  }


  #--------------------------------------------------------------------------------
  # FIND THE MAX OVER THE GRID
  #--------------------------------------------------------------------------------

  if(opt$LocalMax) {

    # Local max.
    if(Grid2d) {

      # Implement local max:
      loc_max_out <- find_local_max(as.array(like))
      indexB <- loc_max_out$row
      indexD <- loc_max_out$col

      # Store value of local max.
      local_max <- list(b = bGrid[indexB],
                      d = dGrid[indexD],
                      like = loc_max_out$value)

      # Record temporary value of maximal likelihood value.
      # It may be replaced to avoid the identification problem.
      max_like <- max(like, na.rm = TRUE)

    } else {

      # Implement local max:
      # Need to pad 1-D array with boundary columns.
      loc_max_out <- find_local_max(as.array(cbind(like - 1,
                                                   like,
                                                   like - 1)))
      indexB <- loc_max_out$row
      indexD <- 1

      # Store value of local max.
      local_max <- list(b = bGrid[indexB],
                        d = rep(NA, length(bGrid[indexB])),
                        like = loc_max_out$value)

      # Record temporary value of maximal likelihood value.
      # It may be replaced to avoid the identification problem.
      max_like <- max(like, na.rm = TRUE)
    }


    # If there is no local max, return global max.
    if(is.null(indexB) | is.null(indexD)) {
      warning('Failure to find local maximum. Returning global maximum instead.',
              'Check whether global maximum is on a boundary and, if so, ',
              'consider expanding the constraints on d and b. ')

      max_like <- max(like, na.rm = TRUE)
      indexB_and_D <- which(like == max_like, arr.ind = TRUE)

      indexB <- indexB_and_D[, 1, drop = FALSE]
      indexD <- indexB_and_D[, 2, drop = FALSE]
    }


  } else {
    # Global max.
    max_like <- max(like, na.rm = TRUE)
    indexB_and_D <- which(like == max(like, na.rm = TRUE), arr.ind = TRUE)

    indexB <- indexB_and_D[, 1, drop = FALSE]
    indexD <- indexB_and_D[, 2, drop = FALSE]
  }


  # Choose value among maxima, if not unique.
  if(length(indexD) > 1 | length(indexB) > 1) {
    # If maximum is not unique, take the index pair corresponding to
    # the highest value of b.
    if(Grid2d) {
      # Sort in ascending order according to indexB.
      # In a 2-D grid, the ordering is unambiguous:
      # sorting in increasing indexB is also increasing in b.
      indBindx <- order(indexB)
      indexB <- indexB[indBindx]
      indexD <- indexD[indBindx]
    } else {
      # Sort in ascending order by b.
      # indexB <- indexB[order(indexB)]

      # Not so fast!
      # In a 1-D grid, the ordering is also unambiguous
      # WHEN THERE ARE NO RESTRICTIONS, aside from d = b:
      # sorting in increasing indexB is also increasing in b.
      # However, when a linear restriction is placed on d and b,
      # the grid is on the 1-D parameter phi,
      # and the relationship with b is either positive or negative,
      # depending on the second element of H in H*[d, b]' + h.
      if(!is.null(opt$R_psi)) {
        sign_H2 <- sign(H[2])
      } else {
        sign_H2 <- 1
      }
      # Then, sort in the appropriate order.
      indexB <- indexB[order(sign_H2*indexB)]

    }

    # Take last value, which is the local maximum with highest value of b.
    indexB <- indexB[length(indexB)]
    indexD <- indexD[length(indexD)] # Works even in 1-D grid, since both ordered by b.

    # Choosing this local maximum avoids the identification problem.
    # max_like <- local_max$like[which.max(local_max$b)]
    # Again, the order depends on the sign of H[2], if it is a restricted model.
    max_like <- local_max$like[which.max(sign_H2*local_max$b)]

    # cat(sprintf('\nWarning, grid search did not find a unique maximum of the log-likelihood function.\n')  )
    message('Grid search did not find a unique maximum of the log-likelihood function.',
            '\nInspect the plot of the log-likelihood function to verify your choice of the LocalMax option.')
  }

  # Translate to d,b.
  if(!is.null(opt$R_psi)) {
    # Translate the parameter estimates.
    dbHatStar <- t(H %*% bGrid[indexB] + h)

    # Translate the grid points to units of the original parameters.
    # This is useful for plotting the likelihood function.
    bGrid_orig <- 0*bGrid
    dGrid_orig <- 0*dGrid
    for (i in 1:length(bGrid)) {
      dbHatStar_i <- t(H %*% bGrid[i] + h)
      dGrid_orig[i] <- dbHatStar_i[1]
      bGrid_orig[i] <- dbHatStar_i[2]
    }

    # Translate local maxima.
    for (i in 1:length(local_max$b)) {
      # Recall that, in this case, bGrid is actually the grid for phi,
      # the 1-D search parameter in the restricted definition of [d, b].
      local_max_db <- matrix(c(local_max$b[i]), nrow = length(local_max$b[i]), ncol = 1)

      dbHatStar_i <- t(H %*% local_max_db + h)
      # These were switched:
      # local_max$b[i] <- dbHatStar_i[1]
      # local_max$d[i] <- dbHatStar_i[2]
      # d is the first parameter in [d, b]; b is the second.
      local_max$d[i] <- dbHatStar_i[1]
      local_max$b[i] <- dbHatStar_i[2]
    }

  } else {
    # No translation necessary.
    dbHatStar <- cbind(dGrid[indexD], bGrid[indexB])
    bGrid_orig <- bGrid
    dGrid_orig <- dGrid
  }



  # Print final progress report.
  if (opt$progress != 0) {
    if ( iterCount == totIters) {
      if (opt$progress == 1) {
        utils::setTxtProgressBar(progbar, value = iterCount)
        message(sprintf('\nProgress : %5.1f%%, b = %4.2f, d = %4.2f, like = %g.',
                        (iterCount/totIters)*100,
                        dbHatStar[2], dbHatStar[1], max(like, na.rm = TRUE) ))
      } else if (opt$progress == 2) {
        message(sprintf('Progress : %5.1f%%, b = %4.2f, d = %4.2f, like = %g.',
                        (iterCount/totIters)*100,
                        dbHatStar[2], dbHatStar[1], max(like, na.rm = TRUE) ))
      }

    }

  }


  #--------------------------------------------------------------------------------
  # STORE THE PARAMETER VALUES
  #--------------------------------------------------------------------------------

  # Output a list of values, including the previous vector \code{params}.
  likeGrid_params <- list(
    params = NA,
    dbHatStar = dbHatStar,
    muHatStar = NA,
    # max_like = max(like, na.rm = TRUE),
    max_like = max_like,
    local_max = NA,
    Grid2d = Grid2d,
    dGrid = dGrid,
    bGrid = bGrid,
    dGrid_orig = dGrid_orig,
    bGrid_orig = bGrid_orig,
    like = like
  )

  # Record local maxima, if required.
  if(opt$LocalMax) {
    likeGrid_params$local_max = local_max
  }


  # Set parameter values for output.
  params <- dbHatStar
  likeGrid_params$params <- params

  # Add level parameter corresponding to max likelihood.
  if(opt$levelParam) {

    muHatStar  <- mu[, indexB, indexD]
    likeGrid_params$muHatStar <- muHatStar


    # Concatenate parameter values.
    params <- matrix(c(params, muHatStar),
                     nrow = 1, ncol = length(params) + length(muHatStar))


    likeGrid_params$params <- params

  }

  # Append additional parameters and set class of output.
  likeGrid_params$k <- k
  likeGrid_params$r <- r
  likeGrid_params$opt <- opt
  class(likeGrid_params) <- 'FCVAR_grid'


  #--------------------------------------------------------------------------------
  # PLOT THE LIKELIHOOD
  #--------------------------------------------------------------------------------


  if (opt$plotLike) {

    graphics::plot(x = likeGrid_params)

  }


  return(likeGrid_params)
}


#' Plot the Likelihood Function for the FCVAR Model
#'
#' \code{plot.FCVAR_grid} plots the likelihood function from \code{FCVARlikeGrid}.
#' \code{FCVARlikeGrid} performs a grid-search optimization
#' by calculating the likelihood function
#' on a grid of candidate parameter values.
#' This function evaluates the likelihood over a grid of values
#' 	for \code{c(d,b)} (or \code{phi}, when  there are constraints on \code{c(d,b)).
#' 	It can be used when parameter estimates are sensitive to
#' 	starting values to give an approximation of the global max which can
#' 	then be used as the starting value in the numerical optimization in
#' 	\code{FCVARestn}.
#'
#' @param x An S3 object of type \code{FCVAR_grid} output from \code{FCVARlikeGrid}.
#' @param y An argument for generic method \code{plot} that is not used in \code{plot.FCVAR_grid}.
#' @param ... Arguments to be passed to methods, such as graphical parameters
#' for the generic plot function.
#' @return NULL
#' @examples
#' opt <- FCVARoptions()
#' opt$gridSearch   <- 0 # Disable grid search in optimization.
#' opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
#' opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
#' opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
#' opt$progress <- 2 # Show progress report on each value of b.
#' # newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
#' likeGrid_params <- FCVARlikeGrid(x, k = 2, r = 1, opt)
#' \dontrun{plot(likeGrid_params)}
#' @note Calls \code{graphics::persp} when \code{x$Grid2d == TRUE} and
#' calls \code{graphics::plot} when \code{x$Grid2d == FALSE}.
#' @family FCVAR auxilliary functions
#' @seealso \code{FCVARoptions} to set default estimation options.
#' \code{plot.FCVAR_grid} plots the likelihood function from \code{FCVARlikeGrid}.
#' @export
#'
plot.FCVAR_grid <- function(x, y = NULL, ...) {

  # Extract parameters.
  Grid2d <- x$Grid2d
  like <- x$like
  dGrid <- x$dGrid
  bGrid <- x$bGrid
  dGrid_orig <- x$dGrid_orig
  bGrid_orig <- x$bGrid_orig
  opt <- x$opt

  # Additional parameters.
  dots <- list(...)


  # if (!is.null(main)) {
  if ('main' %in% names(dots)) {
    main <- dots$main
  } else {
      main <- c('Log-likelihood Function ',
                sprintf('Rank: %d, Lags: %d', x$r, x$k))
  }


  # Plot likelihood depending on dimension of search.
  if(Grid2d) {
    # 2-dimensional plot.

    colors <- grDevices::rainbow(100)
    like2D <- t(like)
    like.facet.center <- (like2D[-1, -1] + like2D[-1, -ncol(like2D)] + like2D[-nrow(like2D), -1] + like2D[-nrow(like2D), -ncol(like2D)])/4
    # Range of the facet center on a 100-scale (number of colors)
    like.facet.range <- cut(like.facet.center, 100)

    # Reduce the size of margins.
    # graphics::par()$mar
    # 5.1 4.1 4.1 2.1
    # bottom, left, top and right margins respectively.
    graphics::par(mar = c(1.1, 1.1, 2.1, 1.1))
    # Reset after:
    # graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

    graphics::persp(dGrid_orig, bGrid_orig,
          like2D,
          xlab = 'd',
          ylab = 'b',
          zlab = '',
          main = main,
          phi = 45, theta = -55,
          r = 1.5,
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
      plot_grid <- bGrid_orig
    } else {
      like_x_label <- 'phi'
      plot_grid <- bGrid
    }

    graphics::plot(plot_grid, like,
         main = main,
         ylab = 'log-likelihood',
         xlab = like_x_label,
         type = 'l',
         col = 'blue',
         lwd = 3)

  }


}


# Find local optima
#
# \code{find_local_max} finds the row and column numbers corresponding
# to the local maxima of a 1- or 2-dimensional array.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param x A 1- or 2-dimensional numerical array of the candidate values
# of the objective function.
# @return A list object with three elements:
# \describe{
#   \item{\code{row}}{An integer row number of the local maxima. }
#   \item{\code{col}}{An integer column number of the local maxima. }
#   \item{\code{value}}{The numeric values of the local maxima. }
# }
#
#
find_local_max <- function(x) {

  # Determine size of input.
  nrows <- dim(x)[1]
  ncols <- dim(x)[2]

  # Remove corner cases.
  if (nrows == 1 & ncols == 1) {
    return(list(row = 1, col = 1, value = x))
  } else if (nrows == 1) {
    x <- t(x)
    nrows <- ncols
    ncols <- 1
  }

  # Initialize matrix of local max indicators.
  loc_max_ind <- matrix(TRUE, nrow = nrows, ncol = ncols)

  # Proceed by ruling out locally dominated points.

  # Check above.
  x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
  x_check[ - 1, ] <- x[ - nrows, ]
  loc_max_ind <- loc_max_ind & (x > x_check)

  # Check below.
  x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
  x_check[ - nrows, ] <- x[ - 1, ]
  loc_max_ind <- loc_max_ind & (x > x_check)


  # Check sides and diagonals only if 2-dimensional matrix.
  if (ncols > 1) {

    # Check left.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ , - 1] <- x[ , - ncols]
    loc_max_ind <- loc_max_ind & (x > x_check)

    # Check right.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ , - ncols] <- x[ , - 1]
    loc_max_ind <- loc_max_ind & (x > x_check)


    # Check top left.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - 1, - 1] <- x[ - nrows, - ncols]
    loc_max_ind <- loc_max_ind & (x > x_check)

    # Check top right.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - 1, - ncols] <- x[ - nrows, - 1]
    loc_max_ind <- loc_max_ind & (x > x_check)


    # Check bottom left.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - nrows, - 1] <- x[ - 1, - ncols]
    loc_max_ind <- loc_max_ind & (x > x_check)

    # Check bottom right.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - nrows, - ncols] <- x[ - 1, - 1]
    loc_max_ind <- loc_max_ind & (x > x_check)

  }


  # Assemble output into a list.
  loc_max_out <- list(
    row = matrix(rep(seq(nrows), ncols),
                 nrow = nrows,
                 ncol = ncols)[loc_max_ind],
    col = matrix(rep(seq(ncols), nrows),
                 nrow = nrows,
                 ncol = ncols, byrow = TRUE)[loc_max_ind],
    value = x[loc_max_ind])

  return(loc_max_out)
}



# Likelihood Function for the Unconstrained FCVAR Model
#
# \code{FCVARlikeMu} calculates the likelihood for the unconstrained FCVAR model
# for a given set of parameter values. It is used by the \code{FCVARlikeGrid}
# function to numerically optimize over the level parameter for given values of
# the fractional parameters.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param y A matrix of variables to be included in the system.
# @param db The orders of fractional integration.
# @param mu The level parameter.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A number \code{like}, the log-likelihood evaluated at the
# specified parameter values.
# @examples
# opt <- FCVARoptions()
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# like <- FCVARlikeMu(mu = colMeans(x), y = x, db = c(1, 1), k = 2, r = 1, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# The \code{FCVARlikeGrid} calls this function to perform a grid search over the
# parameter values.
# @export
#
FCVARlikeMu <- function(mu, y, db, k, r, opt) {


  t <- nrow(y)
  x <- y - matrix(1, nrow = t, ncol = 1) %*% mu


  # Obtain concentrated parameter estimates.
  estimates <- GetParams(x, k, r, db, opt)

  # Calculate value of likelihood function.
  cap_T <- t - opt$N
  p <- ncol(x)
  like <- - cap_T*p/2*( log(2*pi) + 1)  - cap_T/2*log(det(estimates$OmegaHat))

  return(like)
}


# Likelihood Function for the Constrained FCVAR Model
#
# \code{FCVARlike} calculates the likelihood for the constrained FCVAR model
# for a given set of parameter values. This function adjusts the variables with the level parameter,
# if present, and returns the log-likelihood given \code{d} and \code{b}.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param params A vector of parameters \code{d} and \code{b}
# (and \code{mu} if option selected).
# @param x A matrix of variables to be included in the system.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A number \code{like}, the log-likelihood evaluated at the
# specified parameter values.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# FCVARlike(c(results$coeffs$db, results$coeffs$muHat), x, k = 2, r = 1, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# @export
#
FCVARlike <- function(params, x, k, r, opt) {

  # global estimatesTEMP

  cap_T <- nrow(x)
  p <- ncol(x)


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

    y <- x - matrix(1, nrow = cap_T, ncol = p) %*% diag(mu)
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
  # Base this on raw inputs:
  cap_T <- nrow(x) - opt$N
  p <- ncol(x)
  like <- - cap_T*p/2*( log(2*pi) + 1)  - cap_T/2*log(det(estimates$OmegaHat))


  return(like)
}


# Obtain concentrated parameter estimates for the Constrained FCVAR Model
#
# \code{GetEstimates} calculates the concentrated parameter estimates for
# the constrained FCVAR model for a given set of parameter values.
# It is used after estimating the parameters that are numerically optimized
# to obtain the corresponding concentrated parameters.
# Like \code{FCVARlike} , this function adjusts the variables with the level parameter,
# if present, and returns the log-likelihood given \code{d} and \code{b}.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param params A vector of parameters \code{d} and \code{b}
# (and \code{mu} if option selected).
# @param x A matrix of variables to be included in the system.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A list object \code{estimates} containing the estimates,
# including the following parameters:
# \describe{
#   \item{\code{db}}{The orders of fractional integration (taken directly from the input).}
#   \item{\code{alphaHat}}{A \eqn{p x r} matrix of adjustment parameters.}
#   \item{\code{betaHat}}{A \eqn{p x r} matrix of cointegrating vectors.
#       The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.}
#   \item{\code{rhoHat}}{A \eqn{p x 1} vector of restricted constatnts.}
#   \item{\code{piHat}}{A \eqn{p x p} matrix \eqn{\Pi = \alpha \beta'} of long-run levels.}
#   \item{\code{OmegaHat}}{A \eqn{p x p} covariance matrix of the error terms.}
#   \item{\code{GammaHat}}{A ( \eqn{p x kp} matrix \code{cbind(GammaHat1,...,GammaHatk)})
#   of autoregressive coefficients. }
#   \item{\code{muHat}}{A vector of the optimal \code{mu} if level parameter is selected. }
# }
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# GetEstimates(c(results$coeffs$db, results$coeffs$muHat), x, k = 2, r = 1, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARlike} performs the same calculations to obtain the value
# of the likelihood function.
# @export
#
GetEstimates <- function(params, x, k, r, opt) {

  cap_T <- nrow(x)
  p <- ncol(x)


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

    y <- x - matrix(1, nrow = cap_T, ncol = p) %*% diag(mu)
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

  # If level parameter is present, they are the last p parameters in the
  #   params vector
  if (opt$levelParam) {
    estimates$muHat <- mu
  } else {
    estimates$muHat <- NULL
  }

  return(estimates)

}


# Likelihood Function for the FCVAR Model
#
# \code{FCVARlikeFull} calculates the likelihood for the constrained FCVAR model
# for a given set of parameter values.
# This function takes the full set of coefficients from a call
# to \code{FCVARestn} and returns the log-likelihood given \code{d} and \code{b}.
#
# @param x A matrix of variables to be included in the system.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param coeffs A list of coefficients of the FCVAR model.
# It is an element of the list \code{results} returned by \code{FCVARestn},
# without the parameters \code{betaHat} and \code{rhoHat}.
# @param betaHat A \eqn{p x r} matrix of cointegrating vectors.
# The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.
# @param rhoHat A \eqn{p x 1} vector of restricted constatnts.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A number \code{like}, the log-likelihood evaluated at the
# specified parameter values.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# FCVARlikeFull(x, k = 2, r = 1, coeffs = results$coeffs,
#               beta = results$coeffs$betaHat, rho = results$coeffs$rhoHat, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} for the estimation of coefficients in \code{coeffs}.
# @export
#
FCVARlikeFull <- function(x, k, r, coeffs, betaHat, rhoHat, opt) {


  # Add betaHat and rhoHat to the coefficients to get residuals because they
  #   are not used in the Hessian calculation and are missing from the
  #   structure coeffs
  coeffs$betaHat <- betaHat
  coeffs$rhoHat <- rhoHat

  # Obtain residuals.
  epsilon <- GetResiduals(x, k, r, coeffs, opt)


  # Calculate value of likelihood function.
  cap_T <- nrow(x) - opt$N
  p <- ncol(x)
  OmegaHat <- t(epsilon) %*% epsilon/cap_T
  like <- - cap_T*p/2*( log(2*pi) + 1)  - cap_T/2*log(det(OmegaHat))

  return(like)
}


# Transform Data for Regression
#
# \code{TransformData} transforms the raw data by fractional differencing.
# The output is in the format required for regression and
# reduced rank regression.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param x A matrix of variables to be included in the system.
# @param k The number of lags in the system.
# @param db The orders of fractional integration.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A list object \code{Z_array} containing the transformed data,
# including the following parameters:
# \describe{
#   \item{\code{Z0}}{A matrix of data calculated by fractionally differencing
#   \code{x} at differencing order \code{d}.}
#   \item{\code{Z1}}{The matrix \code{x} (augmented with a vector of ones
#   if the model includes a restricted constant term), which is then lagged and stacked
#   and fractionally differenced (with order \eqn{d-b}).}
#   \item{\code{Z2}}{The matrix \code{x} lagged and stacked
#   and fractionally differenced (with order \code{d}).}
#   \item{\code{Z3}}{A column of ones if model includes an
#   unrestricted constant term, otherwise \code{NULL}.}
# }
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# Z_array <- TransformData(x, k = 2, db = results$coeffs$db, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} calls \code{GetParams}, which calls \code{TransformData}
# to estimate the FCVAR model.
# \code{TransformData} in turn calls \code{FracDiff} and \code{Lbk}
# to perform the transformation.
# @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
# "Likelihood inference for a fractionally cointegrated
# vector autoregressive model," Econometrica 80, 2667-2732.
# @references Johansen, S. (1995). "Likelihood-Based Inference
# in Cointegrated Vector Autoregressive Models,"
# New York: Oxford University Press.
# @export
#
TransformData <- function(x, k, db, opt) {

  # Number of initial values and sample size.
  N <- opt$N
  cap_T <- nrow(x) - N

  # Extract parameters from input.
  d <- db[1]
  b <- db[2]

  # Transform data as required.
  Z0 <- FracDiff(x, d)

  Z1 <- x
  # Add a column with ones if model includes a restricted constant term.
  if(opt$rConstant) {
    Z1 <- cbind(x, matrix(1, nrow = N + cap_T, ncol = 1))
  }

  Z1 <- FracDiff(  Lbk( Z1 , b, 1)  ,  d - b )

  Z2 <- FracDiff(  Lbk( x , b, k)  , d)

  # Z3 contains the unrestricted deterministics
  Z3 <- NULL
  # # Add a column with ones if model includes a unrestricted constant term.
  if(opt$unrConstant) {
    Z3 <- matrix(1, nrow = cap_T, ncol = 1)
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


# Calculate Residuals for the FCVAR Model
#
# \code{GetResiduals} calculates residuals for the FCVAR model
# from given parameter values.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param x A matrix of variables to be included in the system.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param coeffs A list of coefficients of the FCVAR model.
# It is an element of the list \code{results} returned by \code{FCVARestn}.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A matrix \code{epsilon} of residuals from FCVAR model estimation
# calculated with the parameter estimates specified in \code{coeffs}.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# epsilon <- GetResiduals(x, k = 2, r = 1, coeffs = results$coeffs, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} to estimate the FCVAR model.
# @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
# "Likelihood inference for a fractionally cointegrated
# vector autoregressive model," Econometrica 80, 2667-2732.
# @export
#
GetResiduals <- function(x, k, r, coeffs, opt) {


  # If level parameter is included, the data must be shifted before
  #   calculating the residuals:
  if (opt$levelParam) {
    cap_T <- nrow(x)
    y <- x - matrix(1, nrow = cap_T, ncol = 1) %*% coeffs$muHat
  } else {
    y <- x
  }

  #--------------------------------------------------------------------------------
  # Transform data
  #--------------------------------------------------------------------------------
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


# Calculate Lag Polynomial in the Fractional Lag Operator
#
# \code{Lbk} calculates a lag polynomial in the fractional lag operator.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param x A matrix of variables to be included in the system.
# @param b The order of fractional differencing.
# @param k The number of lags in the system.
# @return A matrix \code{Lbkx} of the form \eqn{[ Lb^1 x, Lb^2 x, ..., Lb^k x]}
# where \eqn{Lb = 1 - (1-L)^b}.
# The output matrix has the same number of rows as \code{x}
# but \code{k} times as many columns.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# Lbkx <- Lbk(x, b = results$coeffs$db[2], k = 2)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} calls \code{GetParams}, which calls \code{TransformData}
# to estimate the FCVAR model.
# \code{TransformData} in turn calls \code{FracDiff} and \code{Lbk}
# to perform the transformation.
# @export
#
Lbk <- function(x, b, k) {

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
#' In particular, the difference between the out put of the function calls
#' \code{FCVAR::FracDiff(x - mean(x), d = 0.5)}
#' and \code{fracdiff::diffseries(x, d = 0.5)} is numerically small.
#' @export
#'
FracDiff <- function(x, d) {


  if(is.null(x)) {
    dx <- NULL
  } else {


    p <- ncol(x)
    cap_T <- nrow(x)
    np2 <- stats::nextn(2*cap_T - 1, 2)

    k <- 1:(cap_T-1)
    b <- c(1, cumprod((k - d - 1)/k))


    dx <- matrix(0, nrow = cap_T, ncol = p)

    for (i in 1:p) {


      dxi <- stats::fft(stats::fft(c(b, rep(0, np2 - cap_T))) *
                          stats::fft(c(x[, i], rep(0, np2 - cap_T))), inverse = cap_T) / np2

      dx[, i] <- Re(dxi[1:cap_T])

    }



  }

  return(dx)
}


# Count the Number of Free Parameters
#
# \code{GetFreeParams} counts the number of free parameters based on
# 	the number of coefficients to estimate minus the total number of
# 	restrictions. When both \code{alpha} and \code{beta} are restricted,
# 	the rank condition is used to count the free parameters in those two variables.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param p The number of variables in the system.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @param rankJ The rank of a conditioning matrix, as described in
# Boswijk & Doornik (2004, p.447), which is only used if there are
# restrictions imposed on \code{alpha} or \code{beta}, otherwise \code{NULL}.
# @return The number of free parameters \code{fp}.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = NULL)
#
# opt <- FCVARoptions()
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# opt$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = 4)
#
# opt <- FCVARoptions()
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# opt$R_Alpha <- matrix(c(0, 1, 0), nrow = 1, ncol = 3)
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# GetFreeParams(k = 2, r = 1, p = 3, opt = newOpt, rankJ = 4)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn}, \code{HypoTest} and \code{LagSelect} to estimate the FCVAR model
# and use this in the calculation of the degrees of freedom
# for a variety of statistics.
# @references Boswijk, H. P. and J. A. Doornik (2004).
# "Identifying, estimating and testing restricted cointegrated systems:
# An overview," Statistica Neerlandica 58, 440-465.
# @export
#
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


# Calculate the Hessian Matrix
#
# \code{FCVARhess} calculates the Hessian matrix of the
# 	log-likelihood by taking numerical derivatives.
# 	It is used to calculate the standard errors of parameter estimates.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param x A matrix of variables to be included in the system.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param coeffs A list of coefficients of the FCVAR model.
# It is an element of the list \code{results} returned by \code{FCVARestn}.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return The \code{hessian} matrix  of second derivatives of the FCVAR
# log-likelihood function, calculated with the parameter estimates
# specified in \code{coeffs}.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# hessian <- FCVARhess(x, k = 2, r = 1, coeffs = results$coeffs, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} to estimate the FCVAR model and calculate
# standard errors of the estimates.
# @export
#
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


# Collect Parameters into a Vector
#
# \code{SEmat2vecU} transforms the model parameters in matrix
# 	form into a vector.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param coeffs A list of coefficients of the FCVAR model.
# It is an element of the list \code{results} returned by \code{FCVARestn}.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param p The number of variables in the system.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A vector \code{param} of parameters in the FCVAR model.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# params <- SEmat2vecU(coeffs = results$coeffs, k = 2, r = 1, p = 3, opt)
# coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
#
# params <- matrix(seq(25))
# coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
# params <- SEmat2vecU(coeffs = coeffs, k = 2, r = 1, p = 3, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} to estimate the FCVAR model and calculate
# standard errors of the estimates.
# \code{SEmat2vecU} is called by \code{FCVARhess} to sort the parameters
# into a vector to calculate the Hessian matrix.
# \code{SEvec2matU} is a near inverse of \code{SEmat2vecU},
# in the sense that \code{SEvec2matU} obtains only a
# subset of the parameters in \code{results$coeffs}.
#
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


# Extract Parameters from a Vector
#
# \code{SEvec2matU} transforms the vectorized model parameters
# 	into matrices.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param param A vector of parameters in the FCVAR model.
# @param k The number of lags in the system.
# @param r The cointegrating rank.
# @param p The number of variables in the system.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return \code{coeffs}, a list of coefficients of the FCVAR model.
# It has the same form as an element of the list \code{results}
# returned by \code{FCVARestn}.
#
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
# results <- FCVARestn(x, k = 2, r = 1, opt)
# params <- SEmat2vecU(coeffs = results$coeffs, k = 2, r = 1, p = 3, opt)
# coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
#
# params <- matrix(seq(25))
# coeffs <- SEvec2matU(param = params, k = 2, r = 1, p = 3, opt )
# params <- SEmat2vecU(coeffs = coeffs, k = 2, r = 1, p = 3, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} to estimate the FCVAR model and calculate
# standard errors of the estimates.
# \code{SEmat2vecU} is called by \code{FCVARhess} to convert the parameters
# from a vector into the coefficients after calculating the Hessian matrix.
# \code{SEmat2vecU} is a near inverse of \code{SEvec2matU},
# in the sense that \code{SEvec2matU} obtains only a
# subset of the parameters in \code{results$coeffs}.
#
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




# Calculate Restricted Estimates for the FCVAR Model
#
# \code{GetRestrictedParams} calculates restricted estimates of
# cointegration parameters and error variance in the FCVAR model.
# To calculate the optimum, it uses the switching algorithm
# of Boswijk and Doornik (2004, page 455) to optimize over free parameters
# \eqn{\psi} and \eqn{\phi} directly, combined with the line search proposed by
# Doornik (2016, working paper). We translate between  \eqn{(\psi, \phi)} and
# \eqn{(\alpha, \beta)} using the relation of \eqn{R_\alpha vec(\alpha) = 0} and
# \eqn{A\psi = vec(\alpha')}, and \eqn{R_\beta vec(\beta) = r_\beta} and
# \eqn{H\phi + h = vec(\beta)}. Note the transposes.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param beta0 The unrestricted estimate of \code{beta},
# a \eqn{p x r} matrix of cointegrating vectors
# returned from \code{FCVARestn} or \code{GetParams}.
# @param S00 A matrix of product moments,
# calculated from the output of \code{TransformData} in \code{GetParams}.
# @param S01 A matrix of product moments,
# calculated from the output of \code{TransformData} in \code{GetParams}.
# @param S11 A matrix of product moments,
# calculated from the output of \code{TransformData} in \code{GetParams}.
# @param cap_T The number of observations in the sample.
# @param p The number of variables in the system.
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
#
# @return A list object \code{switched_mats} containing the restricted estimates,
# including the following parameters:
# \describe{
#   \item{\code{betaStar}}{A \eqn{p x r} matrix of cointegrating vectors.
#       The \eqn{r x 1} vector \eqn{\beta x_t} is the stationary cointegration relations.}
#   \item{\code{alphaHat}}{A \eqn{p x r} matrix of adjustment parameters.}
#   \item{\code{OmegaHat}}{A \eqn{p x p} covariance matrix of the error terms.}
# }
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# opt$R_Beta <- matrix(c(1, 0, 0), nrow = 1, ncol = 3)
# opt <- FCVARoptionUpdates(opt, p = 3, r = 1) # Need to update restriction matrices.
# betaStar <- matrix(c(-0.95335616, -0.07345676, -0.29277318), nrow = 3)
# S00 <- matrix(c(0.0302086527,  0.001308664,  0.0008200923,
#                 0.0013086640,  0.821417610, -0.1104617893,
#                 0.0008200923, -0.110461789,  0.0272861128), nrow = 3)
# S01 <- matrix(c(-0.0047314320, -0.04488533,  0.006336798,
#                 0.0026708007,  0.17463884, -0.069006455,
#                 -0.0003414163, -0.07110324,  0.022830494), nrow = 3, byrow = TRUE)
# S11 <- matrix(c( 0.061355941,  -0.4109969,  -0.007468716,
#                  -0.410996895,  70.6110313, -15.865097810,
#                  -0.007468716, -15.8650978,   3.992435799), nrow = 3)
# switched_mats <- GetRestrictedParams(betaStar, S00, S01, S11, cap_T = 316, p = 3, opt)
# @family FCVAR auxilliary functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} calls \code{GetParams} to estimate the FCVAR model,
# which in turn calls \code{GetRestrictedParams} if there are restrictions
# imposed on \code{alpha} or \code{beta}.
# @references Boswijk, H. P. and J. A. Doornik (2004).
# "Identifying, estimating and testing restricted cointegrated systems:
# An overview," Statistica Neerlandica 58, 440-465.
# @references Doornik, J. A. (2018).
# "Accelerated estimation of switching algorithms: the cointegrated
# VAR model and other applications,"
# Forthcoming in Scandinavian Journal of Statistics.
# @export
#
GetRestrictedParams <- function(beta0, S00, S01, S11, cap_T, p, opt) {


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
  # iters   <- opt$UncFminOptions$MaxFunEvals
  # Tol     <- opt$UncFminOptions$TolFun
  iters   <- opt$unc_optim_control$maxit
  Tol     <- opt$unc_optim_control$reltol
  conv    <- 0
  i       <- 0

  # Tolerance for entering the line search epsilon_s in Doornik's paper.
  TolSearch <- 0.01

  # Line search parameters.
  lambda <- c(1, 1.2, 2, 4, 8)
  nS <- length(lambda)
  likeSearch <- matrix(NA, nrow = nS, ncol = 1)
  OmegaSearch <- array(0, dim = c(p, p,nS))


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






