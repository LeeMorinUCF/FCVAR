
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


#' Grid Search over Likelihood Values
#'
#' \code{LikeGrid} calculates the likelihood function
#' on a grid of candidate parameter values.
#' This function evaluates the likelihood over a grid of values
#' 	for \code{c(d,b)} (or \code{phi}).
#' 	It can be used when parameter estimates are sensitive to
#' 	starting values to give an approximation of the global max which can
#' 	then be used as the starting value in the numerical optimization in
#' 	\code{FCVARestn}.
#'
#' @param x A matrix of variables to be included in the system.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{EstOptions()}.
#' @return A vector \code{params} of \code{d} and \code{b}
#' (and \code{mu} if level parameter is selected)
#' corresponding to a maximum over the grid of \code{c(d,b)} or \code{phi}.
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' params <- LikeGrid(x, k = 2, r = 1, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{EstOptions} to set default estimation options.
#' @note If \code{opt$LocalMax == 0}, \code{LikeGrid} returns the parameter values
#'       corresponding to the global maximum of the likelihood on the grid.
#'       If \code{opt$LocalMax == 1}, \code{LikeGrid} returns the parameter values for the
#'       local maximum corresponding to the highest value of \code{b}. This
#'       alleviates the identification problem mentioned in Johansen and
#'       Nielsen (2010, section 2.3).
#' @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2010).
#' "Likelihood inference for a nonstationary fractional
#' autoregressive model," Journal of Econometrics 158, 51-66.
#'
LikeGrid <- function(x, k, r, opt) {



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
    dbStep <- 0.02
    # dbStep <- 0.2
  } else {
    Grid2d <- 0
    dbStep <- 0.01
  }



  # If equality restrictions are imposed, need to construct a grid over
  #   phi and obtain db <- H*phi + h.
  if(!is.null(opt$R_psi)) {
    # This set of variables is defined for easy translation between
    #   phi (unrestricted parameter) and d,b (restricted parameters).
    R <- opt$R_psi
    s <- opt$r_psi
    H <- null(R)
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
      cat(sprintf('Now estimating for iteration %d of %d: b = %f.\n ',
                  iterCount, totIters, b))
    }

    # If d>=b (constrained) then search only over values of d that are
    #   greater than or equal to b. Also, perform this operation only if
    #   searching over both d and b.
    if(opt$constrained & Grid2d) {
      # Start at the index that corresponds to the first value in the
      # grid for d that is >= b
      # dStart <-  find(dGrid>=b, 1)
      dStart <- dGrid[dGrid >= b][1]
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

        min_out <- optim(StartVal,
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


      if (opt$progress != 0) {
        # if (toc(lastTic) > opt$updateTime | iterCount == totIters) {
        if ( iterCount == totIters) {
          if (opt$progress == 1) {
            SimNotes <- sprintf('Model: k=%g, r=%g\nb=%4.2f, d=%4.2f, like=%8.2f',
                                k, r, db[2],db[1], like[iB,iD] )
            # waitbar(iterCount/totIters,M_status_bar, SimNotes)
          } else {
            cat(sprintf('Progress : %5.1f%%, b=%4.2f, d=%4.2f, like=%g\n',
                        (iterCount/totIters)*100, db[2], db[1], like[iB,iD] ))
          }

          # lastTic <- tic()
        }

      }


    }


  }

  # # Save the workspace at this point.
  # save.image(file = 'LikeGridData1.RData')
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


  if(length(indexD)>1 | length(indexB)>1) {
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
    cat(sprintf('\nWarning, grid search did not find a unique maximum of the log-likelihood function.\n')  )
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


  params <- dbHatStar

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


    # print('cbind(params, muHatStar) = ')
    # print(cbind(params, muHatStar))
    #
    # params <- cbind(params, muHatStar)

    # Don't ask:
    params <- matrix(c(params, muHatStar),
                     nrow = 1, ncol = length(params) + length(muHatStar))
    # Vectors and matrices and arrays. Oh my!


    # print('params = ')
    # print(params)

  }


  #--------------------------------------------------------------------------------
  # PLOT THE LIKELIHOOD
  #--------------------------------------------------------------------------------


  if (opt$plotLike) {

    # cat(sprintf("Sorry, but I don't feel like plotting the likelihood function right now."))
    # Oh, alright, I'll plot it anyway.

    if(Grid2d) {
      # 2-dimensional plot.

      # Color palette (100 colors)
      col.pal <- colorRampPalette(c("blue", "red"))
      colors <- col.pal(100)
      # height of facets
      like.facet.center <- (like[-1, -1] + like[-1, -ncol(like)] + like[-nrow(like), -1] + like[-nrow(like), -ncol(like)])/4
      # Range of the facet center on a 100-scale (number of colors)
      like.facet.range <- cut(like.facet.center, 100)


      persp(dGrid, bGrid,
            like,
            phi = 45, theta = 45,
            xlab = 'd',
            ylab = 'b',
            main = c('Log-likelihood Function ',
                     sprintf('Rank: %d, Lags: %d', r, k)),
            # col = 'red',
            col = colors[like.facet.range]
      )

    } else {
      # 1-dimensional plot.
      if (is.null(opt$R_psi)) {
        like_x_label <- 'd=b'
      } else {
        like_x_label <- 'phi'
      }

      plot(bGrid, like,
           main = c('Log-likelihood Function ',
                    sprintf('Rank: %d, Lags: %d', r, k)),
           ylab = 'log-likelihood',
           xlab = like_x_label,
           type = 'l',
           col = 'blue',
           lwd = 3)

    }

  }



  return(params)
}


#' Likelihood Function for the Unconstrained FCVAR Model
#'
#' \code{FCVARlikeMu} calculates the likelihood for the unconstrained FCVAR model
#' for a given set of parameter values. It is used by the \code{LikeGrid}
#' function to numerically optimize over the level parameter for given values of
#' the fractional parameters.
#'
#' @param y A matrix of variables to be included in the system.
#' @param db The orders of fractional integration.
#' @param mu The level parameter.
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{EstOptions()}.
#' @return A number \code{like}, the log-likelihood evaluated at the
#' specified parameter values.
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' like <- FCVARlikeMu(y = x, db = c(1, 1), mu = mean(x), k = 2, r = 1, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{EstOptions} to set default estimation options.
#' The \code{LikeGrid} calls this function to perform a grid search over the
#' parameter values.
#'
FCVARlikeMu <- function(mu, y, db, k, r, opt) {


  t <- nrow(y)
  x <- y - matrix(1, nrow = t, ncol = 1) %*% mu

  # print(summary(x))

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
#' @param x A matrix of variables to be included in the system.
#' @param params A vector of parameters \code{d} and \code{b}
#' (and \code{mu} if option selected).
#' @param k The number of lags in the system.
#' @param r The cointegrating rank.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{EstOptions()}.
#' @return A number \code{like}, the log-likelihood evaluated at the
#' specified parameter values.
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' like <- FCVARlike(x, params = c(1, 1), , k = 2, r = 1, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{EstOptions} to set default estimation options.
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
  if(nrow(opt$R_psi) == 1) {
    H <- null(opt$R_psi)
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
  estimatesTEMP <- estimates
  # If level parameter is present, they are the last p parameters in the
  #   params vector
  if (opt$levelParam) {
    estimatesTEMP$muHat <- mu
  } else {
    estimatesTEMP$muHat <- NULL
  }


  # Calculate value of likelihood function.
  T <- nrow(y) - opt$N
  p <- ncol(y)
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates$OmegaHat))

  # Assign the value for the global variable estimatesTEMP.
  # Might adjust this later but it will work for now.
  # print('Search list in FCVARlike():')
  # print(search())
  assign("estimatesTEMP", estimatesTEMP, envir = .GlobalEnv)

  return(like)
}


#' Likelihood Function for the FCVAR Model
#'
#' \code{FullFCVARlike} calculates the likelihood for the constrained FCVAR model
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
#' The \eqn{r x 1} vector \eqn{\betax_t} is the stationary cointegration relations.
#' @param rhoHat A \eqn{p x 1} vector of restricted constatnts.
#' @param opt A list object that stores the chosen estimation options,
#' generated from \code{EstOptions()}.
#' @return A number \code{like}, the log-likelihood evaluated at the
#' specified parameter values.
#' @examples
#' opt <- EstOptions()
#' x <- data(JNP2014)
#' results <- FCVARestn(x, k = 2, r = 1, opt)
#' like <- FullFCVARlike(x, k = 2, r = 1, coeffs = results$coeffs,
#' beta = coeffs$betaHat, rho = coeffs$rhoHat, opt)
#' @family FCVAR auxilliary functions
#' @seealso \code{EstOptions} to set default estimation options.
#' \code{FCVARestn} for the estimation of coefficients in \code{coeffs}.
#'
FullFCVARlike <- function(x, k, r, coeffs, beta, rho, opt) {


  # Add betaHat and rhoHat to the coefficients to get residuals because they
  #   are not used in the Hessian calculation and are missing from the
  #   structure coeffs
  coeffs$betaHat <- beta
  coeffs$rhoHat <- rho

  # Obtain residuals.
  epsilon <- GetResiduals(x, k, r, coeffs, opt)


  # Calculate value of likelihood function.
  T <- nrow(x) - opt$N
  p <- ncol(x)
  OmegaHat <- t(epsilon) %*% epsilon/T
  like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(OmegaHat))

  return(like)
}










