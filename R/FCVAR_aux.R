
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

