################################################################################
# 
# Estimation Options for FCVAR
# 
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 10, 2020
# 
################################################################################
# 
# Written by Michal Popiel and Morten Nielsen (This version 04.11.2016)
#
# DESCRIPTION: This class defines the estimation options used in the FCVAR
# 	estimation procedure and the related programs. Assigning this class
# 	to a variable stores the default properties defined below in that
# 	variable. In addition to the properties, the methods section includes the
# 	function updateRestrictions which performs several checks on the
# 	user-specified options prior to estimation.
# 
################################################################################


################################################################################
# Load packages
################################################################################

# Load pracma package of tools for linear algebra. 
# install.packages('pracma')
library(pracma)
# Dependency: null function. 


################################################################################
# Define list object to contain estimation options.
################################################################################

obj <- list(
  
  #--------------------------------------------------------------------------------
  # ESTIMATION OPTIONS 
  #--------------------------------------------------------------------------------
  
  # Estimation options for unconstrained optimization.
  UncFminOptions = list(MaxFunEvals = 1000, TolFun = 1e-8, 
                         TolX = 1e-8, Display = 'off'),
  
  # Estimation options for constrained optimization.
  ConFminOptions = list(MaxFunEvals = 1000, TolFun = 1e-8, 
                         TolX = 1e-8, Display = 'off', 
                         Algorithm = 'interior-point'),
  
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
  # dbMin<= b <= d <= dbMax. If this option is set to 0 then
  # dbMin<= d <= dbMax and dbMin <= b <= dbMax are imposed separately.
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
  # specified if constrained=0. Furthermore, the grid search is not
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


################################################################################
# Define function to set and test estimation options.
################################################################################
# 
# function newObj = updateRestrictions(obj, p, r)
# Description: This function is called prior to estimation and
# 	performs several checks to verify the input in the estimation
# 	options:
# 	- Check that only one option for deterministics is chosen,
# 	- Check that starting values have been specified correctly,
# 	- Update dbMin, dbMax, based on user-specified options,
# 	- Check for appropriate dimensions and redundancies in the
# 		restriction matrices R_psi, R_Alpha and R_Beta.
# Input = obj (estimation option class variable)
#         p (number of variables in the system)
#         r (cointegrating rank)
# Output = newObj (updated object)
# Dependencies: called by FCVARestn().
# 
################################################################################


updateRestrictions <- function(obj, p, r) {
  
  #--------------------------------------------------------------------------------
  # Deterministics
  #--------------------------------------------------------------------------------
  
  # Check for redundant deterministic specifications.
  if(obj$levelParam && obj$rConstant) {
    
    obj$rConstant <- 0
    print('\nWarning, both restricted constant and level parameter selected.\n')
    print('Only level parameter will be included.\n')
    
  }
  
  #--------------------------------------------------------------------------------
  # Fractional Parameters d,b
  #--------------------------------------------------------------------------------
  
  # Adjust starting values.
  
  # If d=b is imposed but the starting values are different, then
  # adjust them so that they are the same for faster and more
  # accurate computation.
  if(obj$restrictDB & length(obj$db0) > 1 &
     obj$db0[1] != obj$db0[2]) {
    
    obj$db0 <- c(obj$db0[1], obj$db0[1])
    print('\nWarning, d=b imposed but starting values are different for d,b.\n')
    print(sprintf('db0 set to (%g, %g).\n', obj$db0[1], obj$db0[2]))
    
  }
  
  # if only one starting value is specified, then d0 = b0.
  if(length(obj$db0) < 2) {
    obj$db0 <- c(obj$db0, obj$db0)
  }
  
  
  # Check if too many parameters specified in dbMin/dbMax
  if ( length(obj$dbMin) > 2 | length(obj$dbMax) > 2 ) {
    
    
    print('\nWARNING: Too many parameters specified in dbMin or dbMax.\n')
    print('\nOnly the first two elements will be used to set bounds in estimation.\n')
    # Cut off extra bounds
    if ( length(obj$dbMin) > 2 ) {
      obj$dbMin <- obj$dbMin[1:2]
    }
    if ( length(obj$dbMax) > 2 ) {
      obj$dbMax <- obj$dbMax[1:2]
    }
    
  }
  
  
  
  
  # Assign the same min/max values for d,b if only one set of
  # bounds is provided. This is mostly for backwards
  # compatibility with previous versions.
  # Minimum
  if( length(obj$dbMin) == 1 ) {
    obj$dbMin <- c(obj$dbMin, obj$dbMin)
    # Notification to user.
    print('\nWarning: minimum value specified for d only.\n')
    print(sprintf('Min value of b set to %2.4f\n', obj$dbMin[2]))
  }
  else {
    # Set as row vector, regardless of input.
    # obj$dbMin <- reshape(obj$dbMin,1,2)
    # Not needed since R "recycles". 
    obj$dbMin <- matrix(obj$dbMin, nrow = 1, ncol = 2)
  }
  # Maximum
  if(length(obj$dbMax) == 1) {
    obj$dbMax <- c(obj$dbMax, obj$dbMax)
    # Notification to user.
    print('\nWarning: maximum value specified for d only.\n')
    print(sprintf('Max value of b set to %2.4f\n', obj$dbMax[2]))
  }
  else {
    # Set as row vector, regardless of input.
    # obj$dbMax <- reshape(obj$dbMax,1,2)
    # Not needed since R "recycles". 
    obj$dbMax <- matrix(obj$dbMax, nrow = 1, ncol = 2)
  }
  
  
  
  # Check for consistency among min and max values for fractional
  # parameters.
  if( obj$dbMin[1] > obj$dbMax[1] ) {
    print('\nWarning, min > max for bounds on d.\n')
    stop('Invalid bounds inmposed on d.')
  }
  if( obj$dbMin[2] > obj$dbMax[2] ) {
    print('\nWarning, min > max for bounds on b.\n')
    stop('Invalid bounds inmposed on b.')
  }
  
  
  # If both inequality constraints and d>=b is imposed then
  # inconsistencies could arise in the optimization. Furthermore,
  # grid search is not equipped to deal with other types of
  # inequalities.
  if(!is.null(obj$C_db) & obj$constrained) {
    
    print('\nWarning, inequality restrictions imposed and constrained selected.\n')
    print('Restrictions in C_db override constrained options.\n')
    obj$constrained <- 0
    
    if(obj$gridSearch) {
      print('Grid search has been switched off.\n')
      obj$gridSearch <- 0
    }
    
  }
  
  
  # If restrictDB is selected, then only restrictions on d
  # are allowed.
  if (!is.null(obj$R_psi)) {
    
    # First check if restriction is valid
    if(obj$restrictDB & ((obj$R_psi[1,1] == 0) | 
                          (obj$R_psi[1,2] != 0)) ) {
      print('\nError in R_psi. When restrictDB = 1, only ')
      print('restrictions on d can be imposed.\n')
      stop('Invalid restriction for d=b model.')
    }
    
    # If non-zero restrictions haven't been specified by
    # the user, then set r_psi to zero's.
    if(is.null(obj$r_psi)) {
      obj$r_psi <- 0
    }
    
  }
  
  
  # If d=b is imposed, add it to the other equality restrictions on d,b.
  if(obj$restrictDB) {
    
    # Check conformability: 
    # obj$R_psi <- [obj$R_psi [1 -1]]
    # Stacking a matrix on a new row:
    obj$R_psi <- rbind(obj$R_psi, c(1, -1))
    obj$r_psi <- c(obj$r_psi, 0)
    if(obj$constrained) {
      print('\nNote: Redundant options. Both constrained (d>=b) and restrict (d=b) selected.')
      print('\n Only d=b imposed.\n')
      # Turn off (d>=b) constraint.
      obj$constrained <- 0
    }
    
  }
  
  
  
  # Check restrictions on fractional parameters.
  if (is.null(obj$R_psi)) {
    
    # Ensure that parameter space for fractional parameters is non-empty.
    if(obj$dbMax[1] < obj$dbMin[2] & obj$constrained) {
      print('\nWarning: Redefine restrictions on fractional parameters.\n')
      stop('Empty parameter space for (d,b).')
    }
    
  } else {
    # Ensure that parameter space for fractional parameters is non-empty.
    UB_LB_bounds <- GetBounds(obj)
    UB <- UB_LB_bounds$UB
    LB <- UB_LB_bounds$LB
    
    if(LB > UB) {
      print('\nWarning: Redefine restrictions on fractional parameters.\n')
      stop('Empty parameter space for (d,b).')
    }
    
    # Check if grid search is necessary.
    if ( (nrow(obj$R_psi) > 1 & obj$gridSearch) ) {
      print('\nd and b are exactly identified by imposed restrictions\n')
      print('so grid search has been turned off.\n')
      obj$gridSearch = 0
    }
    
    # Check for redundancies.
    if(qr(obj$R_psi)$rank < nrow(obj$R_psi)) {
      print('\nWARNING: R_psi has reduced rank!\n')
      print('\nRedefine R_psi with linearly independent restrictions.\n')
      stop('Redundant restrictions in R_psi matrix')
    }
    
  }
  
  #--------------------------------------------------------------------------------
  # Alpha and Beta
  #--------------------------------------------------------------------------------
  
  # Define p1 to be number of rows of betaStar.
  p1 = p + obj$rConstant
  
  # --- Alpha --- #
  if(!is.null(obj$R_Alpha)) {
    
    # Check if restricted parameters actually exist.
    if(r == 0) {
      print('\nWARNING: Cannot impose restrictions on alpha if r=0!\n')
      print('Either remove the restriction (set R_Alpha = NULL)  or increase rank (set r>0).\n')
      stop('Imposing restrictions on empty parameters')
    }
    
    # Check if column length of R_Alpha matches the number of
    # parameters.
    if(ncol(obj$R_Alpha) != p %*% r) {
      print('\nWARNING: The length of R_Alpha does not match the number of parameters!\n')
      print('Please respecify R_Alpha so that the number of columns is p*r.\n')
      stop('Restriction misspecification')
    }
    
    # Check for redundancies.
    if(qr(obj$R_Alpha)$rank < nrow(obj$R_Alpha)) {
      print('\nWARNING: R_Alpha has reduced rank!\n')
      print('\nRedefine R_Alpha with linearly independent restrictions.\n')
      stop('Redundant restrictions in R_Alpha matrix')
    }
    end
    # Check if user has imposed non-homogeneous alpha restrictions.
    if(!is.null(obj$r_Alpha) & (obj$r_Alpha != 0)) {
      print('\nWARNING: r_Alpha contains non-homogeneous restrictions (r_alpha non-zero).\n')
      print('All alpha restrictions have been made homogeneous.\n')
    }
    
    obj$r_Alpha <- matrix(0, nrow = nrow(obj$R_Alpha), ncol = 1)
    
  }
  
  
  # --- Beta --- #
  if(!is.null(obj$R_Beta)) {
    
    
    # Check if restricted parameters actually exist.
    if(r == 0) {
      print('\nWARNING: Cannot impose restrictions on Beta if r=0!\n')
      print('Either remove the restriction (set R_Beta = NULL)  or increase rank (set r>0).\n')
      stop('Imposing restrictions on empty parameters')
    }
    
    # Check if column length of R_Beta matches the number of
    # parameters (note p1 not p).
    if(ncol(obj$R_Beta) != p1 %*% r) {
      print('\nWARNING: The length of R_Beta does not match the number of parameters!\n')
      print('Please respecify R_Beta so that the number of columns is p1*r.\n')
      stop('Restriction misspecification')
    }
    
    # Check for redundancies.
    if(qr(obj$R_Beta)$rank < nrow(obj$R_Beta)) {
      print('\nWARNING: R_Beta has reduced rank!\n')
      print('\nRedefine R_Beta with linearly independent restrictions.\n')
      stop('Redundant restrictions in R_Beta matrix')
    }
    
    
    # If non-zero restrictions haven't been specified by
    # the user, then set the r_Beta vector to zero's.
    if(is.null(obj$r_Beta)) {
      obj$r_Beta <- matrix(0, nrow = nrow(obj$R_Beta), ncol = 1)
    }
    else {
      # If user has specified restrictions, then check if
      # the dimensions of LHS and RHS match. Note that if
      # a restricted constant is being estimated, then it
      # needs to be accounted for in the restrictions.
      if(nrow(obj$r_Beta) != nrow(obj$R_Beta)) {
        print('\nWARNING: Row dimensions of R_Beta and r_Beta do not match!\n')
        print('Please redefine these matrices so that dimensions match.\n')
        print('Note: if a restricted constant has been included, \n')
        print('the Beta matrix has an additional row.\n')
        stop('All restrictions must be specified in the r_Beta variable')
      }
      
    }
    
    
  }
  
  
  
  # Return the updated object of estimation options.
  newObj <- obj
  return(newObj)
  
}


################################################################################
# Define function to obtain bounds on parameters.
################################################################################
# 
# function [ UB, LB ] <- GetBounds(opt)
# Written by Michal Popiel and Morten Nielsen (This version 04.11.2016)
# 
# DESCRIPTION: This function obtains upper and lower bounds on d,b or on 
#   phi, given by db <- H*phi + h. 
#
# Input  <- opt   (object containing estimation options)
# Output <- UB (a 2x1 or 1x1 upper bound for db or phi)
#          LB (a 2x1 or 1x1 lower bound for db or phi)
# 
################################################################################


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
    H <- null(R)
    h <- t(R) %*% inv(R*t(R)) %*% s 
    
    UB <- NULL
    LB <- NULL
    
    # If there are two restrictions imposed, both d and b are restricted
    # and the upper and lower bounds are set to their restricted values.
    if(nrow(R) == 2) {
      UB <- t(h)
      LB <- t(h)
    }
    end
    
    # If there is 1 restriction, then there is 1 free parameter.
    if(size(R,1) == 1) {
      
      # Calculate endpoints of the grid from dbMin and dbMax to free
      # parameter phi. Note that since the null space can be either
      # positive or negative, the following conditional statements are
      # needed.
      if(h[1]>0) {
        phiMin1 <- (opt$dbMin(1) - h[1]) / h[1]
        phiMax1 <- (opt$dbMax(1) - h[1]) / h[1]
      } else if(h[1]<0) {
        phiMin1 <- (opt$dbMax(1) - h[1]) / h[1]
        phiMax1 <- (opt$dbMin(1) - h[1]) / h[1]
      }
      else {
        phiMin1 <- NA
        phiMax1 <- NA
      }
      end
      if(h[2]>0) {
        phiMin2 <- (opt$dbMin(2) - h[2]) / h[2]
        phiMax2 <- (opt$dbMax(2) - h[2]) / h[2]
      } else if(h[2]<0) {
        phiMin2 <- (opt$dbMax(2) - h[2]) / h[2]
        phiMax2 <- (opt$dbMin(2) - h[2]) / h[2]
      } else {
        phiMin2 <- NA
        phiMax2 <- NA
      }
      end
      
      # Take into account the condition d>=b, if required.
      if(opt$constrained) {
        if (h[1]>h[2]) {
          phiMin3 <- (h[2]-h[1])/(h[1] - h[2])
          phiMax3 <- NA
        }
        else {
          phiMax3 <- (h[2]-h[1])/(h[1] - h[2])
          phiMin3 <- NA
        }
        end
      }
      else {
        phiMin3 <- NA
        phiMax3 <- NA
      }
      end
      
      LB <- max(max(phiMin1, phiMin2), phiMin3)
      UB <- min(min(phiMax1, phiMax2), phiMax3)    
      end
      
      
    }
    
  }
  
  UB_LB_bounds <- list(UB = UB, LB = LB)
  
  return(UB_LB_bounds)
}

#_________________________________________________________________________


################################################################################
# End
################################################################################



################################################################################
# Testing null() function
################################################################################
# 
# 
# install.packages('pracma')
# library(pracma)
# 
# a <- matrix(c(
#   16,     2,     3,    13,
#   5,    11,    10,     8,
#   9,     7,     6,    12,
#   4,    14,    15,     1), 4, 4, byrow = TRUE)
# 
# null(a)
# # [,1]
# # [1,] -0.2236068
# # [2,] -0.6708204
# # [3,]  0.6708204
# # [4,]  0.2236068
# 
# 
# a <- matrix(c(
#   1,8,15,67,
#   7,    14,    16,     3), nrow = 2, ncol = 4, byrow = TRUE)
# 
# 
# null(a)
# # [,1]        [,2]
# # [1,] -0.29767375 -0.89698149
# # [2,] -0.63974275  0.43972988
# # [3,]  0.70442146  0.01565751
# # [4,] -0.07687621 -0.04262269
# 
# # Match the values for the matlab function.
# # https://www.mathworks.com/help/matlab/ref/null.html
# 


