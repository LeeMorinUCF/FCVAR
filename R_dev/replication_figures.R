
# Old functions for plotting FCVAR results for JSS draft.


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
#' \dontrun{plot.FCVARlikeGrid(likeGrid_params, k = 2, r = 1, newOpt, main = 'default')}
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



# Example of function with dots as arguments.
test_dots <- function(x, ...) {

  dots <- list(...)
  print(names(dots))
  if ('main' %in% names(dots)) {
    main <- dots$main
    print('main = ')
    print(main)
  } else {
    print('main not specified.')
  }

  return(100*x^2)
}


test_dots(x = 4)
test_dots(x = 5, main = 'bananas')
test_dots(x = 6, y = 7, fruit = 'oranges')
