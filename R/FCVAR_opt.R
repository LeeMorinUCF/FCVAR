


#' Set Estimation Options
#'
#' \code{FCVARoptions} defines the estimation options used in the FCVAR
#'   estimation procedure and the related programs.
#'
#' @param ... a list of arguments to set to values other than the default settings.
#' See argument names in return value below.
#' @return An S3 object of class \code{FCVAR_opt} that stores the default estimation options,
#' which includes the following parameters:
#' \describe{
#'   \item{\code{unc_optim_control}}{A list of options in the form of the argument \code{control}
#'   in the \code{optim} function for \emph{unconstrained} optimization of the likelihood function
#'   over the fractional integration parameters.
#'   This is also used in the switching algorithm employed when linear constraints are imposed on
#'   the cointegrating relations \code{beta} or the adjustment coefficients \code{alpha},
#'   so it must at least contain the arguments \code{maxit} and \code{reltol},
#'   since it uses those parameters.}
#'   \item{\code{con_optim_control}}{A list of options in the form of the argument \code{control}
#'   in either the \code{optim} or the \code{constrOptim} function for \emph{constrained} optimization
#'   of the likelihood function over the fractional integration parameters,
#'   using the 'L-BFGS-B' algorithm.
#'   It must at least contain the arguments \code{maxit} and \code{pgtol}.}
#'   \item{\code{LineSearch}}{Indicator for conducting a line search optimization within
#'   the switching algorithm when optimizing over constraints on the cointegrating relations \eqn{\beta}
#'   or the adjustment coefficients \eqn{\alpha}. See Doornik (2018, Section 2.2) for details.}
#'   \item{\code{LocalMax}}{Indicator to select the local maximum with the highest value of \code{b}
#'   when there are multiple local optima. This is meant to alleviate the identification problem discussed
#'    in Johansen and Nielsen (2010, Section 2.3)  and Carlini and de Magistris (2019).
#'    When \code{LocalMax <- 0}, the optimization returns the values of \code{d} and \code{b}
#'    corresponding to the global optimum.}
#'   \item{\code{dbMax}}{Upper bound for the fractional integration parameters \code{d}, \code{b}.}
#'   \item{\code{dbMin}}{Lower bound for the fractional integration parameters \code{d}, \code{b}.}
#'   \item{\code{db0}}{The starting values for optimization of the fractional integration parameters \code{d}, \code{b}.}
#'   \item{\code{constrained}}{Indicator to impose restriction \code{dbMax >= d >= b >= dbMin}.}
#'   \item{\code{restrictDB}}{Indicator to impose restriction \code{d = b}.}
#'   \item{\code{N}}{The number of initial values: the observations to condition upon.}
#'   \item{\code{unrConstant}}{Indicator to include an unrestricted constant.}
#'   \item{\code{rConstant}}{Indicator to include an restricted constant.}
#'   \item{\code{levelParam}}{Indicator to include level parameter.}
#'   \item{\code{C_db}}{CHECK whether still used.}
#'   \item{\code{c_db}}{CHECK whether still used.}
#'   \item{\code{UB_db}}{An upper bound on the fractional integration parameters \code{d} and \code{b},
#'   after transforming the parameters to account for any restrictions imposed.}
#'   \item{\code{LB_db}}{An lower bound on the fractional integration parameters \code{d} and \code{b},
#'   after transforming the parameters to account for any restrictions imposed.}
#'   \item{\code{R_psi}}{A matrix for defining restrictions on the fractional integration
#'   parameters \code{d} and \code{b}, of the form \eqn{R_{\psi}(d, b)' = r_{\psi}}.}
#'   \item{\code{r_psi}}{A vector for defining restrictions on the fractional integration
#'   parameters \code{d} and \code{b}, of the form \eqn{R_{\psi}(d, b)' = r_{\psi}}.}
#'   \item{\code{R_Alpha}}{A matrix for defining restrictions on the adjustment coefficients
#'   of the form \eqn{R_{\alpha}\alpha = r_{\alpha}}.}
#'   \item{\code{r_Alpha}}{A vector for defining restrictions on the adjustment coefficients
#'   of the form \eqn{R_{\alpha}\alpha = r_{\alpha}}.}
#'   \item{\code{R_Beta}}{A matrix for defining restrictions on the cointegrating relations
#'   of the form \eqn{R_{\beta}\beta = r_{\beta}}.}
#'   \item{\code{r_Beta}}{A vector for defining restrictions on the cointegrating relations
#'   of the form \eqn{R_{\beta}\beta = r_{\beta}}.}
#'   \item{\code{print2screen}}{Indicator to print output to screen.}
#'   \item{\code{printGammas}}{Indicator to print estimates and standard errors on autoregressive
#'   coefficients \eqn{\Gamma_i, i = i, ..., k}.}
#'   \item{\code{printRoots}}{Indicator to print roots of characteristic polynomial.}
#'   \item{\code{plotRoots}}{Indicator to plot roots of characteristic polynomial.}
#'   \item{\code{CalcSE}}{Indicator to calculate the standard errors. It is used when displaying results.}
#'   \item{\code{hess_delta}}{Size of increment for numerical calculation of derivatives of the likelihood
#'   function for numerical calculation of the Hessian matrix. The default is \code{10^(-4)},
#'   which works well in practice to balance errors between precision and truncation.}
#'   \item{\code{gridSearch}}{Indicator to perform a grid search for the optimization
#'   over the fractional integration parameters, for more accurate estimation.
#'   This will make estimation take longer.}
#'   \item{\code{dbStep1D}}{The step size for the grid search over the fractional integration parameters
#'   for the 1-dimensional grid search (such as when restrictions are imposed between \code{d} and \code{b}.). }
#'   \item{\code{dbStep2D}}{The step size for the grid search over the fractional integration parameters
#'   for the 2-dimensional grid search.}
#'   \item{\code{plotLike}}{Indicator to plot the likelihood (only if \code{gridSearch <- 1}).}
#'   \item{\code{progress}}{Show a waitbar for a progress indicator for the grid search.}
#'   \item{\code{updateTime}}{How often progress is updated in the waitbar for the grid search (in seconds).}
#' }
#' @examples
#' opt <- FCVARoptions()
#' opt <- FCVARoptions(
#'     gridSearch   = 0, # Disable grid search in optimization.
#'     dbMin        = c(0.01, 0.01), # Set lower bound for d,b.
#'     dbMax        = c(2.00, 2.00), # Set upper bound for d,b.
#'     constrained  = 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
#' )
#' @family FCVAR estimation functions
#' @seealso \code{FCVARoptionUpdates} to set and test estimation options for validity and compatibility.
#' \code{FCVARestn} for use of these options in estimation.
#' @references Doornik, J. A. (2018) "Accelerated Estimation of Switching Algorithms:
#' The Cointegrated VAR Model and Other Applications."
#' Scandinavian Journal of Statistics, Volume 45, Issue 2.
#' @references Johansen, S\enc{ø}{o}ren, and Morten \enc{Ø}{O}rregaard Nielsen (2010) "Likelihood inference for a nonstationary
#' fractional autoregressive model." Journal of Econometrics 158, 51–66.
#' @references Carlini, F., and P. S. de Magistris (2019) "On the identification of fractionally cointegrated
#' VAR models with the F(d) condition." Journal of Business & Economic Statistics 37(1), 134–146.
#' @export
#'
FCVARoptions <- function(...) {

  # Take list of arguments passed to FCVARoptions, if any.
  opt_dots <- list(...)
  # Any valid arguments here will overwrite default options.
  # Arguments of fixed size (esp. scalars) will be overwritten at the end,
  # and variable-sized arguments are assigned within the list definition.

  # Create list of default options.
  opt <- list(

    #--------------------------------------------------------------------------------
    # ESTIMATION OPTIONS
    #--------------------------------------------------------------------------------

    # Estimation options for unconstrained optimization.
    # UncFminOptions = list(MaxFunEvals = 1000,
    #                       TolFun = 1e-8 #, # Not all are used in R.
    #                       # TolX = 1e-8,
    #                       # Display = 'off'
    # ),
    # unc_optim_control = list(maxit = 1000,
    #                          reltol = 1e-8 #, # Not all are used in R.
    #                          # TolX = 1e-8,
    #                          # Display = 'off'
    # ),
    unc_optim_control = if ('unc_optim_control' %in% names(opt_dots)) {
      opt_dots$unc_optim_control
    } else {
      list(maxit = 1000,
           reltol = 1e-8)
    },

    # Estimation options for constrained optimization.
    # ConFminOptions = list(MaxFunEvals = 1000,
    #                       TolFun = 1e-8,
    #                       # TolX = 1e-8,
    #                       # Display = 'off',
    #                       # Algorithm = 'interior-point',
    #                       Algorithm = 'L-BFGS-B'),
    # con_optim_control = list(maxit = 1000,
    #                          pgtol = 1e-8),
    con_optim_control = if ('con_optim_control' %in% names(opt_dots)) {
      opt_dots$con_optim_control
    } else {
      list(maxit = 1000,
           pgtol = 1e-8)
    },
    # L-BFGS-B is limited-memory BFGS with bounds.

    # Activate live search for switching algorithm in restricted model
    # estimation.
    LineSearch = 1,

    # Variation on the grid search to find local max.
    LocalMax   = 1,

    # Set upper and lower bound for d,b parameters.
    dbMax = if ('dbMax' %in% names(opt_dots)) {opt_dots$dbMax} else {2},
    dbMin = if ('dbMin' %in% names(opt_dots)) {opt_dots$dbMin} else {0.01},

    # Starting value for optimization.
    db0 = c(1, 1),

    # Constrain parameter space for d,b. If set to 1, then impose
    # dbMin <= b <= d <= dbMax. If this option is set to 0 then
    # dbMin <= d <= dbMax and dbMin <= b <= dbMax are imposed separately.
    constrained = 1,

    # If restrictDB = 1, then d=b is imposed. If = 0, the restriction
    # is not imposed.
    restrictDB = 1,

    # Number of observations reserved for initial values to be
    # conditioned upon in the estimation.
    N = 0,


    #--------------------------------------------------------------------------------
    # DETERMINISTICS
    #--------------------------------------------------------------------------------

    # Unrestricted constant.
    unrConstant = 0,

    # Restricted constant.
    rConstant   = 0,

    # Level parameter.
    levelParam  = 1,


    #--------------------------------------------------------------------------------
    # RESTRICTIONS
    #--------------------------------------------------------------------------------

    # Inequality constraints on parameters d and b.
    # Specified as C_db * [d b] <= c_db
    # Note: these restrictions are non-standard and should only be
    # specified if constrained = 0. Furthermore, the grid search is not
    # equipped to handle these types of restrictions.
    C_db = if ('C_db' %in% names(opt_dots)) {opt_dots$C_db} else {NULL},
    c_db = if ('c_db' %in% names(opt_dots)) {opt_dots$c_db} else {NULL},

    # Upper and lower bounds for parameters d and b.
    # Note: these options are set automatically in the
    # FCVARoptionUpdates function below based on inputs of 'dbMax',
    # and 'dbMin'.
    UB_db = if ('UB_db' %in% names(opt_dots)) {opt_dots$UB_db} else {NULL},
    LB_db = if ('LB_db' %in% names(opt_dots)) {opt_dots$LB_db} else {NULL},

    # Equality constraints on parameters d and b.
    # Specified as R_psi * [d b] = r_psi
    R_psi = if ('R_psi' %in% names(opt_dots)) {opt_dots$R_psi} else {NULL},
    r_psi = if ('r_psi' %in% names(opt_dots)) {opt_dots$r_psi} else {NULL},


    # Restrictions on Alpha matrix.
    # Specified as R_Alpha*vec(Alpha) = r_Alpha
    # Note: r_Alpha can only have 0s
    R_Alpha = if ('R_Alpha' %in% names(opt_dots)) {opt_dots$R_Alpha} else {NULL},
    r_Alpha = if ('r_Alpha' %in% names(opt_dots)) {opt_dots$r_Alpha} else {NULL},

    # Restrictions on Beta matrix.
    # Specified as R_Beta*vec(Beta) = r_Beta
    # Note: r_Beta can have non-zero elements.
    R_Beta = if ('R_Beta' %in% names(opt_dots)) {opt_dots$R_Beta} else {NULL},
    r_Beta = if ('r_Beta' %in% names(opt_dots)) {opt_dots$r_Beta} else {NULL},

    #--------------------------------------------------------------------------------
    # OUTPUT OPTIONS
    #--------------------------------------------------------------------------------

    # If set to 0, print2screen prevents all output from being printed.
    print2screen = 1,

    # If set to 0, printGammas prevents the estimated coefficients in
    # the Gamma matrix from being printed.
    printGammas = 1,

    # If set to 0, printRoots prevents the roots of the
    # characteristic polynomial from being printed.
    printRoots = 1,

    # If set to 0, plotRoots prevents the roots of the
    # characteristic polynomial from being plotted.
    plotRoots = 1,


    #--------------------------------------------------------------------------------
    # Switch for calculating SEs
    #--------------------------------------------------------------------------------

    # Calculate SEs. This option should not be changed by the user.
    #   Calculating standard errors is omitted in the LagSelection.m
    #   RankTests.m functions to speed up estimation.
    CalcSE = 1,

    # Standard errors are calculated numerically, by taking numerical centered
    # derivatives of the likelihood function.
    # Set increment for changes along each parameter value.
    hess_delta = 10^(-4),


    #--------------------------------------------------------------------------------
    # Grid search
    #--------------------------------------------------------------------------------

    # gridSearch set to 1 makes the program evaluate the likelihood
    # over a grid of (d,b) and choose the max from the grid as the
    # starting value for the numerical optimization. It is computationally
    # costly to perform but sometimes yields more accurate results.
    gridSearch = 1,

    # Set the size of the increments in the grid in fractional parameters d and b.
    # If linear constraints are imposed, the 1-dimensional value applies.
    # For a search over both d and b, the 2-dimensional value applies.
    dbStep1D = 0.01,
    dbStep2D = 0.02,

    # plotLike makes the program plot the likelihood function over the
    # grid of (d,b) when the grid search is selected.
    plotLike = 1,

    # progress prints the progress of the grid search in
    # either a waitbar window (=1) or to the commandline (=2) or not at
    # all (=0).
    progress = 1,

    # If progress ~= 0, print progress every updateTime seconds.
    updateTime = 5


  )

  # Replace any default options with arguments passed through opt_dots (...).
  # Skip those already assigned (arguments in matrix form).
  matrix_args <- c('unc_optim_control', 'con_optim_control',
                   'dbMax', 'dbMin', 'C_db', 'c_db', 'R_psi', 'r_psi',
                   'R_Alpha', 'r_Alpha', 'R_Beta', 'r_Beta')
  dots_arg_list <- names(opt_dots)[!(names(opt_dots) %in% matrix_args)]

  for (arg_name in dots_arg_list) {
    opt[arg_name] <- unname(unlist(opt_dots[arg_name]))
  }


  # Define class of object.
  class(opt) <- 'FCVAR_opt'

  return(opt)
}


# Update Estimation Options for FCVAR Restrictions
#
# A function to set and test estimation options for validity and compatibility.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
# This function is called prior to estimation and
# 	performs several checks to verify the input in the estimation
# 	options:
# \itemize{
#   \item Check that only one option for deterministics is chosen,
#   \item Check that starting values have been specified correctly,
#   \item Update \code{dbMin}, \code{dbMax}, based on user-specified options,
#   \item Check for appropriate dimensions and redundancies in the restriction matrices \code{R_psi}, \code{R_Alpha} and \code{R_Beta}.
# }
#
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @param p The number of variables in the system.
# @param r The cointegrating rank.
# @return An S3 object of class \code{FCVAR_opt} that stores the default estimation options.
# It is a revised version of the output of \code{FCVARoptions()}.
# @examples
# opt <- FCVARoptions()
# opt$gridSearch   <- 0 # Disable grid search in optimization.
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
# newOpt <- FCVARoptionUpdates(opt, p = 3, r = 1)
# @family FCVAR estimation functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} calls this function at the start of each estimation to verify
# validity of options.
# @export
#
FCVARoptionUpdates <- function(opt, p, r) {

  #--------------------------------------------------------------------------------
  # Deterministics
  #--------------------------------------------------------------------------------

  # Check for redundant deterministic specifications.
  if(opt$levelParam & opt$rConstant) {

    opt$rConstant <- 0
    warning("Both restricted constant and level parameter selected.\n",
            "Only level parameter will be included.\n",
            "To avoid this warning message, please choose one of the above options.")

  }

  #--------------------------------------------------------------------------------
  # Fractional Parameters d,b
  #--------------------------------------------------------------------------------

  # Adjust starting values.

  # If d=b is imposed but the starting values are different, then
  # adjust them so that they are the same for faster and more
  # accurate computation.
  if(opt$restrictDB & length(opt$db0) > 1 &
     abs(opt$db0[1] - opt$db0[2]) > 10^-6) {

    opt$db0 <- c(opt$db0[1], opt$db0[1])
    warning("Restriction d = b imposed but starting values are different for d and b.\n",
            sprintf('Starting values db0 set to c(%g, %g).\n', opt$db0[1], opt$db0[2]))

  }

  # if only one starting value is specified, then d0 = b0.
  if(length(opt$db0) < 2) {
    opt$db0 <- c(opt$db0, opt$db0)
  }


  # Check if too many parameters specified in dbMin/dbMax
  if ( length(opt$dbMin) > 2 | length(opt$dbMax) > 2 ) {

    warning("Too many parameters specified in dbMin or dbMax.\n",
            "Only the first two elements will be used to set bounds in estimation.")

    # Cut off extra bounds
    if ( length(opt$dbMin) > 2 ) {
      opt$dbMin <- opt$dbMin[1:2]
    }
    if ( length(opt$dbMax) > 2 ) {
      opt$dbMax <- opt$dbMax[1:2]
    }

  }




  # Assign the same min/max values for d,b if only one set of
  # bounds is provided. This is mostly for backwards
  # compatibility with previous versions.
  # Minimum
  if( length(opt$dbMin) == 1 ) {
    opt$dbMin <- c(opt$dbMin, opt$dbMin)
    warning("Minimum value specified for d only.\n",
            sprintf('Minimum value of b set to %2.4f\n', opt$dbMin[2]))
  } else {
    # Set as row vector, regardless of input.
    opt$dbMin <- matrix(opt$dbMin, nrow = 1, ncol = 2)
  }
  # Maximum
  if(length(opt$dbMax) == 1) {
    opt$dbMax <- c(opt$dbMax, opt$dbMax)
    warning("Maximum value specified for d only.\n",
            sprintf('Maximum value of b set to %2.4f.\n', opt$dbMax[2]))
  } else {
    # Set as row vector, regardless of input.
    opt$dbMax <- matrix(opt$dbMax, nrow = 1, ncol = 2)
  }



  # Check for consistency among min and max values for fractional
  # parameters.
  if( opt$dbMin[1] > opt$dbMax[1] ) {
    stop("Invalid bounds imposed on d.\n",
         "Minimum > maximum for bounds on d, leaving parameter space empty.")
  }
  if( opt$dbMin[2] > opt$dbMax[2] ) {
    stop("Invalid bounds imposed on b.\n",
         "Minimum > maximum for bounds on b, leaving parameter space empty.")
  }


  # If both inequality constraints and d>=b is imposed then
  # inconsistencies could arise in the optimization. Furthermore,
  # grid search is not equipped to deal with other types of
  # inequalities.
  if(!is.null(opt$C_db) & opt$constrained) {

    warning("Inequality restrictions imposed and constrained model selected.\n",
            "Restrictions in C_db override constrained options.")
    opt$constrained <- 0

    if(opt$gridSearch) {
      message("Grid search has been switched off.")
      opt$gridSearch <- 0
    }

  }


  # If restrictDB is selected, then only restrictions on d
  # are allowed.
  if (!is.null(opt$R_psi)) {

    # First check if restriction is valid
    if (opt$restrictDB & ((opt$R_psi[1,1] == 0) |
                          (opt$R_psi[1,2] != 0)) ) {
      stop("Invalid restriction for d = b model.\n",
           "When restrictDB = 1, only restrictions on d can be imposed.\n",
           "Adjust option R_psi by removing nonzero entries on b.")
    }

    # If non-zero restrictions haven't been specified by
    # the user, then set r_psi to zeros.
    if(is.null(opt$r_psi)) {
      opt$r_psi <- 0
    }

  }


  # If d=b is imposed, add it to the other equality restrictions on d,b.
  if(opt$restrictDB) {

    # Stacking a matrix on a new row:
    # Stacking because the equality restriction might be an additional
    # restriction on d and b.
    opt$R_psi <- rbind(opt$R_psi, c(1, -1))
    opt$r_psi <- c(opt$r_psi, 0)

    if(opt$constrained) {

      warning("Redundant options constraining d and b.\n",
              "Both constrained (d >= b) and restrict (d = b) selected.\n",
              "Only d = b is imposed.")
      # Turn off (d>=b) constraint.
      opt$constrained <- 0
    }

  }


  # Check restrictions on fractional parameters.
  if (is.null(opt$R_psi)) {

    # Ensure that parameter space for fractional parameters is non-empty.
    if(opt$dbMax[1] < opt$dbMin[2] & opt$constrained) {
      stop("Empty parameter space for (d, b).\n",
           "Redefine restrictions on fractional parameters.")
    }

  } else {
    # Ensure that parameter space for fractional parameters is non-empty.
    UB_LB_bounds <- GetBounds(opt)
    UB <- UB_LB_bounds$UB
    LB <- UB_LB_bounds$LB


    # if (LB > UB) {
    if (any(LB > UB)) {
      # cat(sprintf('\nWarning: Redefine restrictions on fractional parameters.\n'))
      # stop('Empty parameter space for (d,b).')
      stop("Empty parameter space for (d, b).\n",
           "Redefine restrictions on fractional parameters.")
    }

    # Check if grid search is necessary.
    if ( (nrow(opt$R_psi) > 1 & opt$gridSearch) ) {
      warning("Fractional differencing parameters d and b\n",
              "are exactly identified by imposed restrictions\n",
              "so grid search has been turned off.")
      opt$gridSearch = 0
    }

    # Check for redundancies.
    if(qr(opt$R_psi)$rank < nrow(opt$R_psi)) {
      stop("Redundant restrictions in R_psi matrix.\n",
           "R_psi has reduced rank!",
           "Redefine R_psi with linearly independent restrictions.")
    }

  }

  #--------------------------------------------------------------------------------
  # Alpha and Beta
  #--------------------------------------------------------------------------------

  # Define p1 to be number of rows of betaStar.
  p1 = p + opt$rConstant

  # --- Alpha --- #
  if(!is.null(opt$R_Alpha)) {

    # Check if restricted parameters actually exist.
    if(r == 0) {
      stop("Imposing restrictions on empty parameters.\n",
           "Cannot impose restrictions on alpha if r = 0!\n",
           "Either remove the restriction (set R_Alpha = NULL)  or increase rank (set r > 0).")
    }

    # Check if column length of R_Alpha matches the number of
    # parameters.
    if(ncol(opt$R_Alpha) != p * r) {
      stop("Restriction misspecification.\n",
           "The length of R_Alpha does not match the number of parameters!\n",
           "Respecify R_Alpha so that the number of columns is p*r.")
    }

    # Check for redundancies.
    if(qr(opt$R_Alpha)$rank < nrow(opt$R_Alpha)) {
      stop("Redundant restrictions in R_Alpha matrix.\n",
           "R_Alpha has reduced rank!",
           "Redefine R_Alpha with linearly independent restrictions.")
    }

    # Check if user has imposed non-homogeneous alpha restrictions.
    if(!is.null(opt$r_Alpha) && (opt$r_Alpha != 0)) {
      cat(sprintf('\nWARNING: r_Alpha contains non-homogeneous restrictions (r_alpha non-zero).\n'))
      cat(sprintf('All alpha restrictions have been made homogeneous.\n'))
      warning("r_Alpha contains non-homogeneous restrictions (r_alpha non-zero).\n",
              "All alpha restrictions have been made homogeneous.")
    }

    opt$r_Alpha <- matrix(0, nrow = nrow(opt$R_Alpha), ncol = 1)

  }


  # --- Beta --- #
  if(!is.null(opt$R_Beta)) {


    # Check if restricted parameters actually exist.
    if(r == 0) {
      stop("Imposing restrictions on empty parameters.\n",
           "Cannot impose restrictions on beta if r = 0!\n",
           "Either remove the restriction (set R_Beta = NULL)  or increase rank (set r > 0).")
    }

    # Check if column length of R_Beta matches the number of
    # parameters (note p1 not p).
    if(ncol(opt$R_Beta) != p1 * r) {
      stop("Restriction misspecification.\n",
           "The length of R_Beta does not match the number of parameters!\n",
           "Respecify R_Beta so that the number of columns is p1*r.")
    }

    # Check for redundancies.
    if(qr(opt$R_Beta)$rank < nrow(opt$R_Beta)) {
      stop("Redundant restrictions in R_Beta matrix.\n",
           "R_Beta has reduced rank!",
           "Redefine R_Beta with linearly independent restrictions.")
    }


    # If non-zero restrictions haven't been specified by
    # the user, then set the r_Beta vector to zero's.
    if(is.null(opt$r_Beta)) {
      opt$r_Beta <- matrix(0, nrow = nrow(opt$R_Beta), ncol = 1)
    }
    else {
      # If user has specified restrictions, then check if
      # the dimensions of LHS and RHS match. Note that if
      # a restricted constant is being estimated, then it
      # needs to be accounted for in the restrictions.
      if(nrow(opt$r_Beta) != nrow(opt$R_Beta)) {
        stop("Row dimensions of R_Beta and r_Beta do not match!\n",
             "Redefine these matrices so that dimensions match.\n",
             "Note: if a restricted constant has been included, \n",
             "the Beta matrix has an additional row.\n",
             "All restrictions must also specify the r_Beta variable.")
      }

    }


  }


  # Return the updated object of estimation options.
  newOpt <- opt
  return(newOpt)

}


# Obtain Bounds on FCVAR Parameters
#
# \code{GetBounds()} obtains bounds on fractional integration parameters in the FCVAR model.
#   This function obtains upper and lower bounds on \code{d}, \code{b} or on
#   \code{phi}, which is given by \code{db <- H*phi + h},
#   using the restrictions implied by \code{H} and \code{h}.
# Note that Roxygen comments are excluded to keep this function internal.
# However, the contents of Roxygen comments are shown below for those who read the scripts.
#
# @param opt An S3 object of class \code{FCVAR_opt} that stores the chosen estimation options,
# generated from \code{FCVARoptions()}.
# @return A list object \code{UB_LB_bounds} containing the upper (\code{UB}) and lower (\code{LB})  bounds on either \code{db} or \code{phi}.
# These vectors are either 1- or 2-dimensional depending on whether a restriction on \code{d} and \code{b} is imposed or not.
# @examples
# opt <- FCVARoptions()
# UB_LB_bounds <- GetBounds(opt)
#
# opt <- FCVARoptions()
# opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
# opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
# UB_LB_bounds <- GetBounds(opt)
# @family FCVAR estimation functions
# @seealso \code{FCVARoptions} to set default estimation options.
# \code{FCVARestn} calls this function at the start of each estimation to specify any bounds on fractional integration parameters.
# @export
#
GetBounds <- function(opt) {


  if(is.null(opt$R_psi)) {
    # R_psi empty. Upper and lower bounds are the max and min values input
    # by the user. Both different and same bounds are permitted for d, b.
    # Maximum
    UB <- opt$dbMax
    # Minimum
    LB <- opt$dbMin

  } else {



    # This set of variables is defined for easy translation between
    # phi (unrestricted parameter) and d,b (restricted parameters).
    R <- opt$R_psi
    s <- opt$r_psi
    H <- pracma::null(R)

    h <- t(R) %*% solve(R %*% t(R)) %*% s


    UB <- NULL
    LB <- NULL

    # If there are two restrictions imposed, both d and b are restricted
    # and the upper and lower bounds are set to their restricted values.
    if(nrow(R) == 2) {
      UB <- t(h)
      LB <- t(h)
    }


    # If there is 1 restriction, then there is 1 free parameter.
    if(nrow(R) == 1) {


      # Calculate endpoints of the grid from dbMin and dbMax to free
      # parameter phi. Note that since the null space can be either
      # positive or negative, the following conditional statements are
      # needed.
      if(H[1] > 0) {
        phiMin1 <- (opt$dbMin[1] - h[1]) / H[1]
        phiMax1 <- (opt$dbMax[1] - h[1]) / H[1]
      } else if(H[1] < 0) {
        phiMin1 <- (opt$dbMax[1] - h[1]) / H[1]
        phiMax1 <- (opt$dbMin[1] - h[1]) / H[1]
      } else {
        phiMin1 <- NA
        phiMax1 <- NA
      }

      if(H[2] > 0) {
        phiMin2 <- (opt$dbMin[2] - h[2]) / H[2]
        phiMax2 <- (opt$dbMax[2] - h[2]) / H[2]
      } else if(H[2] < 0) {
        phiMin2 <- (opt$dbMax[2] - h[2]) / H[2]
        phiMax2 <- (opt$dbMin[2] - h[2]) / H[2]
      } else {
        phiMin2 <- NA
        phiMax2 <- NA
      }


      # Take into account the condition d>=b, if required.
      if(opt$constrained) {
        if (H[1] > H[2]) {
          phiMin3 <- (h[2] - h[1])/(H[1] - H[2])
          phiMax3 <- NA
        }
        else {
          phiMax3 <- (h[2] - h[1])/(H[1] - H[2])
          phiMin3 <- NA
        }

      }
      else {
        phiMin3 <- NA
        phiMax3 <- NA
      }


      LB <- max(c(phiMin1, phiMin2, phiMin3), na.rm = TRUE)
      UB <- min(c(phiMax1, phiMax2, phiMax3), na.rm = TRUE)



    }

  }

  UB_LB_bounds <- list(UB = UB, LB = LB)

  return(UB_LB_bounds)
}
