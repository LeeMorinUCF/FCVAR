


#' Default Estimation Options
#'
#' \code{EstOptions} defines the estimation options used in the FCVAR
#'   estimation procedure and the related programs.
#'
#' @return A list object \code{opt} that stores the default estimation options.
#' @examples
#' EstOptions()
#' @family FCVAR option functions
#' @seealso \code{updateRestrictions} to set and test estimation options for validity and compatibility.
#' @export
#'
#'
EstOptions <- function() {

  opt <- list(

    #--------------------------------------------------------------------------------
    # ESTIMATION OPTIONS
    #--------------------------------------------------------------------------------

    # Estimation options for unconstrained optimization.
    UncFminOptions = list(MaxFunEvals = 1000,
                          TolFun = 1e-8 #, # Not all are used in R.
                          # TolX = 1e-8,
                          # Display = 'off'
    ),

    # Estimation options for constrained optimization.
    ConFminOptions = list(MaxFunEvals = 1000,
                          TolFun = 1e-8,
                          # TolX = 1e-8,
                          # Display = 'off',
                          # Algorithm = 'interior-point',
                          Algorithm = 'L-BFGS-B'),
    # L-BFGS-B is limited-memory BFGS with bounds.

    # Activate live search for switching algorithm in restricted model
    # estimation.
    LineSearch = 1,

    # Variation on the grid search to find local max.
    LocalMax   = 1,

    # Set upper and lower bound for d,b parameters.
    dbMax = 2,
    dbMin = 0.01,

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
    C_db = NULL,
    c_db = NULL,

    # Upper and lower bounds for parameters d and b.
    # Note: these options are set automatically in the
    # updateRestrictions function below based on inputs of 'dbMax',
    # and 'dbMin'.
    UB_db = NULL,
    LB_db = NULL,

    # Equality constraints on parameters d and b.
    # Specified as R_psi * [d b] = r_psi
    R_psi = NULL,
    r_psi = NULL,


    # Restrictions on Alpha matrix.
    # Specified as R_Alpha*vec(Alpha) = r_Alpha
    # Note: r_Alpha can only have 0s
    R_Alpha = NULL,
    r_Alpha = NULL,

    # Restrictions on Beta matrix.
    # Specified as R_Beta*vec(Beta) = r_Beta
    # Note: r_Beta can have non-zero elements.
    R_Beta = NULL,
    r_Beta = NULL,

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
    # Grid search
    #--------------------------------------------------------------------------------

    # gridSearch set to 1 makes the program evaluate the likelihood
    # over a grid of (d,b) and choose the max from the grid as the
    # starting value for the numerical optimization. It is computationally
    # costly to perform but sometimes yields more accurate results.
    gridSearch = 1,

    # plotLike makes the program plot the likelihood function over the
    # grid of (d,b) when the grid search is selected.
    plotLike = 1,

    # progress prints the progress of the grid search in
    # either a waitbar window (=1) or to the commandline (=2) or not at
    # all (=0).
    progress = 1,

    # If progress ~= 0, print progress every updateTime seconds.
    updateTime = 5,

    #--------------------------------------------------------------------------------
    # Set location path with program name
    #--------------------------------------------------------------------------------

    # Location path with program name of fracdist program, if installed,
    #   for calculation of P-values, see MacKinnon and Nielsen (2014, JAE).
    # If this is not installed, the option can be left as-is or blank.
    # Note: use both single (outside) and double quotes (inside). This
    #   is especially important if path name has spaces.

    # Linux example:
    progLoc = '"/usr/bin/fdpval"',

    # Windows example:
    #   program located in folder fdpval in current directory
    # progLoc = '".\fdpval\fdpval"'


    #--------------------------------------------------------------------------------
    # Switch for calculating SEs
    #--------------------------------------------------------------------------------

    # Calculate SEs. This option should not be changed by the user.
    #   Calculating standard errors is omitted in the LagSelection.m
    #   RankTests.m functions to speed up estimation.
    CalcSE = 1



  )


  return(opt)
}


