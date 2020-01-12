################################################################################
# 
# Main Estimation Function for FCVAR
# 
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 11, 2020
# 
################################################################################
# 
# Written by Michal Popiel and Morten Nielsen (This version 01.10.2020)
#
# DESCRIPTION: The main FCVAR estimation function. 
# 
################################################################################


################################################################################
# Load packages
################################################################################

# Load packages. 
# install.packages('name_of_package')
# library(name_of_package)



################################################################################
# Define functions.
################################################################################



################################################################################
# Define function to estimate FCVAR model
################################################################################
# 
# function [ results ] <- FCVARestn(x,k,r,opt)
# Written by Michal Popiel and Morten Nielsen (This version 04.09.2016)
# 
# DESCRIPTION: This function performs estimation of the FCVAR system. It is
# 	the main function in the program with several nested functions, each
# 	described below. It estimates the model parameters, calculates the
# 	standard errors and the number of free parameters, obtains the residuals
# 	and the roots of the characteristic polynomial, and prints the output.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- results (a Matlab structure containing estimation results)
#            - results$startVals     (Starting values used for optimization)
#            - results$options       (Estimation options)
#            - results$like          (Model log-likelihood)
#            - results$coeffs        (Parameter estimates)
#            - results$rankJ         (Rank of Jacobian for
#					identification condition)
#            - results$fp            (Number of free parameters)
#            - results$SE            (Standard errors)
#            - results$NegInvHessian (Negative of inverse Hessian matrix)
#            - results$Residuals     (Model residuals)
#            - results$cPolyRoots    (Roots of characteristic polynomial)
# 
################################################################################

FCVARestn <- function(x,k,r,opt) {
  
  # global estimatesTEMP # Different in R.
  
  # --- Preliminary steps --- %
  T <- nrow(x) - opt$N # number of observations
  p <- ncol(x)         # number of variables
  
  # Update options based on initial user input.
  opt <- updateRestrictions(opt, p, r)
  
  # Clear previous instances of coefficient estimates. This is done
  #   because estimatesTemp is a global structure that will not be cleared
  #   automatically if estimation is interrupted.
  estimatesTEMP <- NULL
  
  
  #--------------------------------------------------------------------------------
  # GRID SEARCH
  #--------------------------------------------------------------------------------
  
  # Hide all warnings for grid search.
  # warning('off', 'all')
  
  # Perform grid search and store results as starting values for
  #   numerical optimization below.
  if(opt$gridSearch) {
    print(sprintf('\nRunning grid search over likelihood for k=%g, r=%g.\n', 
            k,r))
    print(sprintf('This computation can be slow.\n'))
    print(sprintf('Set opt$gridSearch <- 0 to skip it.\n'))
    opt$db0 <- LikeGrid(x, k, r, opt)
    # Change upper and lower bounds to limit the search to some small
    #   interval around the starting values.
    opt$UB_db[1:2] <- pmin(opt$db0[1:2] + c(0.1, 0.1), opt$dbMax)
    opt$LB_db[1:2] <- pmax(opt$db0[1:2] - c(0.1, 0.1), opt$dbMin)
  }
  else {
    # Call to GetBounds returns upper/lower bounds for (d,b) or
    #  depending on whether or not restrictions have been imposed.
    UB_LB_bounds <- GetBounds(obj)
    opt$UB_db <- UB_LB_bounds$UB
    opt$LB_db <- UB_LB_bounds$LB
    
  }
  
  
  
  # Turn warnings back on for main estimation.
  # warning('on','all')
  
  # Hide warnings related to algorithm changes in numerical optimization.
  # warning('off','optim:fminunc:SwitchingMethod')
  # warning('off','optimL:fminunc:SwitchingMethod')
  
  # Hide warnings related to singular matrix in numerical optimization.
  # warning('off','MATLAB:nearlySingularMatrix')
  
  
  
  #--------------------------------------------------------------------------------
  # ESTIMATION
  #--------------------------------------------------------------------------------
  
  
  # Store equality restrictions, inequality restrictions, upper and lower
  #   bounds, and starting values for estimation because they need to be
  #   adjusted based on the presence of level parameter or d=b restriction.
  Rpsi <- opt$R_psi
  rpsi <- opt$r_psi
  startVals <- opt$db0[1:2] # starting values for mu get added later.
  Cdb <- opt$C_db
  cdb <- opt$c_db
  
  # If Rpsi is empty, then optimization is over (d,b), otherwise it is
  #  over phi. If it is over phi, need to make adjustments to startVals
  #  and make Cdb and cdb empty.
  
  
  
  if(size(Rpsi,1)==1) {
    
    H_psi <- null(Rpsi)
    
    # Need to back out phi from given db0
    startVals <- solve(t(H_psi) %*% H_psi) %*% t(H_psi) %*% t(startVals[1:2])
    
    if(opt$gridSearch) {
      # Translate from d,b to phi.
      UB <- solve(t(H_psi) %*% H_psi) %*% t(H_psi) %*% t(opt$UB_db)
      LB <- solve(t(H_psi) %*% H_psi) %*% t(H_psi) %*% t(opt$LB_db)
    }
    else {
      # Otherwise GetBounds returns the values in terms of phi.
      UB <- opt$UB_db
      LB <- opt$LB_db
    }
    
    # Add warning about turning these off if non-empty?
    Cdb <- NULL
    cdb <- NULL
    
    # Turn off equality restrictions if only one restriction is imposed
    #  because they will be imposed inside the profile likelihood
    #  function that is being maximized.
    
    Rpsi <- NULL
    rpsi <- NULL
    
    
  } else {
    
    UB <- opt$UB_db
    LB <- opt$LB_db
    
  }
  
  
  
  
  # If estimation involves level parameter, make appropriate
  #   adjustments so that the dimensions of the restrictions correspond to
  #   the number of parameters being estimated, i.e. fill with zeros.
  if(opt$levelParam) {
    
    
    # If there are equality restrictions add coefficients of 0 to the
    #  level parameter.
    if(!is.null(Rpsi)) {
      # Check conformability:
      Rpsi <- rbind(Rpsi, matrix(0, nrow = nrow(Rpsi), ncol = p))
    }
    
    # If there are inequality restrictions add coefficients of 0 to the
    #  level parameter.
    if(!is.null(Cdb)) {
      # Check conformability:
      Cdb <- rbind(Cdb, matrix(0, nrow = nrow(Cdb), ncol = p))
    }
    
    # If grid search has not been used to set all initial values, add p
    #   starting values, set as the first observations on each variable,
    #   to the startVals variable to account for the level parameter.
    if(opt$gridSearch == 0) {
      # Check conformability:
      startVals <- c(startVals, x[1, ])
    }
    else {
      # Check conformability:
      startVals <- c(startVals, opt$db0[3:length(opt$db0)])
    }
    
    
    # Level parameter is unrestricted, but the length of UB and LB
    #   needs to match the number of parameters in optimization.
    UB <- c(UB, matrix(1,nrow = 1, ncol = p)*Inf)
    LB <- c(LB, -matrix(1,nrow = 1, ncol = p)*Inf)
    
    
  }
  
  
  
  if(size(opt$R_psi)==2) {
    
    
    # d,b are exactly identified by the linear restrictions and Rpsi is
    #  invertible. We use opt$R_psi here because Rpsi is adjusted
    #  depending on the presence of level parameters. Transpose is
    #  necessary to match input/output of other cases.
    dbTemp <- t(t(opt$R_psi) %*% solve(opt$R_psi*t(opt$R_psi)) %*% opt$r_psi)
    y <- x
    
    if(opt$levelParam) {
      # Optimize over level parameter for given (d,b).
      StartVal <- y[1, ]
      
      # [ muHat, maxLike, ! ] 
      # <- fminunc(@( params ) -FCVARlikeMu(x, dbTemp, params, k, r, opt), 
      #            StartVal, opt$UncFminOptions )
      
      # Need to implement optimization correctly. 
      min_out <- optim_unc(-FCVARlikeMu(params, x, dbTemp, k, r, opt), 
                           startVal, opt$UncFminOptions)
      
      muHat <- min_out$par
      maxLike <- min_out$f
      
      
      y <- x - matrix(1, nrow = T+opt$N, ncol = p) %*% diag(muHat)
    }
    else {
      maxLike <- -FCVARlike(dbTemp, y, k, r, opt)
    }
    
    
    # Obtain concentrated parameter estimates.
    estimates <- GetParams(y, k, r, dbTemp, opt)
    
    # Storing the estimates in a global structure here allows us to skip a
    #   call to GetParams after optimization to recover the coefficient
    #   estimates
    estimatesTEMP <- estimates
    # If level parameter is present, they are the last p parameters in the
    #   params vector
    if (opt$levelParam) {
      estimatesTEMP$muHat <- muHat
    }
    else {
      estimatesTEMP$muHat <- NULL
    }
    
    
  } else {
    # [ !, maxLike, ! ] 
    # <- fmincon(@( params ) -FCVARlike(x, params, k, r, opt), 
    #            startVals, Cdb, cdb, Rpsi, rpsi, LB, UB, [], opt$ConFminOptions )
    
    # Need to implement optimization correctly. 
    min_out <- optim_con(-FCVARlike(params, x, k, r, opt), 
                     startVals, Cdb, cdb, Rpsi, rpsi, LB, UB, opt$ConFminOptions)
    
    maxLike <- min_out$f
  }
  

  
  
  
  
  #--------------------------------------------------------------------------------
  # Store outputs
  #--------------------------------------------------------------------------------
  
  # Store the updated estimation options. 
  results$startVals <- startVals
  results$options <- opt
  results$options$UB_db <- UB
  results$options$LB_db <- LB
  
  
  # Adjust the sign of the likelihood and store the results
  maxLike <- -maxLike
  results$like <- maxLike
  
  # Coefficients are taken from a global defined in the likelihood
  #   function
  results$coeffs <- estimatesTEMP
  
  
  #--------------------------------------------------------------------------------
  # CHECK RANK CONDITION
  #--------------------------------------------------------------------------------
  
  
  
  p1 <- p + opt$rConstant
  rankJ <- NULL # initialize the rank
  
  if(r > 0) {
    # If rank is zero then Alpha and Beta are empty
    
    if(is.null(opt$R_Beta)) {
      H_beta <- diag(p1*r)
    } else {
      H_beta <- null(opt$R_Beta)
    }
    
    
    # We use the commutation matrix K_pr to transform vec(A) into vec(A'), # '
    #   see Magnus & Neudecker (1988, p. 47, eqn (1)).
    Ip <- diag(p)
    
    # Need to make sure this builds in the proper order: 
    Kpr <- matrix(kronecker(Ip, diag(r)), nrow = p*r, ncol = p*r)
    if(is.null(opt$R_Alpha)) {
      A <- Kpr %*% diag(p*r)
    }
    else {
      A <- null(opt$R_Alpha %*% solve(Kpr))
    }
    
    
    rA <- ncol(A) # number of free parameters in alpha
    rH <- ncol(H_beta) # number of free parameters in beta (including constant)
    
    # Following Boswijk & Doornik (2004, p.447) identification condition
    # Check conformability: 
    kronA <- kronecker(eye(p), 
                       rbind(results$coeffs$betaHat, 
                             results$coeffs$rhoHat)) %*% rbind(A, zeros(p*r,rH))
    kronH <- kronecker(results$coeffs$alphaHat, 
                       diag(p1)) %*% rbind(matrix(0, nrow = p1*r, ncol = rA), H_beta)
    rankJ <- qr(kronA + kronH)$rank
    
    results$rankJ <- rankJ
  } 
  
  
  #--------------------------------------------------------------------------------
  # CHECK RANKS OF ALPHA AND BETA
  #--------------------------------------------------------------------------------
  
  
  # Check that alpha and beta have full rank to ensure that restrictions
  #   do not reduce their rank.
  if(qr(results$coeffs$alphaHat)$rank < r) {
    print(sprintf('\nWarning: Alpha hat has rank less than r!\n'))
  }
  
  
  if( qr(results$coeffs$betaHat)$rank < r) {
    print(sprintf('\nWarning: Beta hat has rank less than r!\n'))
  }
  
  
  #--------------------------------------------------------------------------------
  # FREE PARAMETERS
  #--------------------------------------------------------------------------------
  
  # Compute the number of free parameters in addition to those in alpha
  #   and beta.
  fp <- FreeParams(k, r, p, opt, rankJ)
  # Store the result.
  results$fp <- fp
  
  
  #--------------------------------------------------------------------------------
  # STANDARD ERRORS
  #--------------------------------------------------------------------------------
  
  
  if(opt$CalcSE) {
    
    # If any restrictions have been imposed, the Hessian matrix must be
    #   adjusted to account for them.
    
    if( !is.null(opt$R_Alpha) | !is.null(opt$R_psi) ) {
      
      # Create R matrix with all restrictions.
      
      # Count the number of restrictions on d,b. Note: opt$R_psi already
      #  contains restrict DB, so the size() is only reliable if
      #  it's turned off. 
      if(!is.null(opt$R_psi)) {
        if(opt$restrictDB)
          rowDB <- nrow(opt$R_psi) - 1
        else
          rowDB <- nrow(opt$R_psi)
        end
      } else {
        # Otherwise d,b are unrestricted.
        rowDB <- 0
      }
        end
        
        # Number of restrictions on alpha.
        rowA  <- nrow(opt$R_Alpha)
        
        # Count the variables.
        colDB <- 1 + !opt$restrictDB
        colA <- p*r
        colG <- p*p*k
        colMu <- opt$levelParam*p
        colRh <- opt$unrConstant*p
        # Length of vec(estimated coefficients in Hessian).
        R_cols  <- colDB + colMu + colRh + colA + colG
        
        # The restriction matrix will have rows equal to the number of
        #   restrictions.
        R_rows <- rowDB + rowA
        
        R <- matrix(0, nrow = R_rows, ncol = R_cols)
        
        # Fill in the matrix R.
        
        # Start with restrictions on (d,b) and note that if the model
        #   with d=b is being estimated, only the first column of R_psi is
        #   considered. In that case, if there are zeros found in the first
        #   column, the user is asked to rewrite the restriction.
        if(rowDB > 0) {
          if(opt$restrictDB) {
            # If the model d=b is being estimated, only one
            #  restriction can be imposed and that is on d.
            R[1:rowDB,1:colDB] <- opt$R_psi[1,1]
          }
          else {
            R[1:rowDB,1:colDB] <- opt$R_psi
          }
          end
        }
        end
        # Put the R_Alpha matrix into the appropriate place in R.
        if(!is.null(opt$R_Alpha)) {
          R[1+rowDB: rowDB + rowA, 
            1 + colDB + colMu + colRh: 
              colDB + colMu + colRh + colA] <- opt$R_Alpha
        }
        end
        
        # Calculate unrestricted Hessian.
        H <- FCVARhess(x, k, r, results$coeffs, opt)
        
        
        # Calculate the restricted Hessian.
        Q <- -solve(H) + solve(H) %*% t(R) %*% 
          solve(R %*% solve(H) %*% t(R)) %*% R %*% solve(H)
      


    } else {
      # Model is unrestricted.
      H <- FCVARhess(x, k, r, results$coeffs, opt)
      Q <- -solve(H)
    }
    end
    
    
  } else {
    NumCoeffs <- length(SEmat2vecU(results$coeffs, k, r, p, opt))
    Q <- matrix(0, nrow = NumCoeffs, ncol = NumCoeffs)
  }
  end
  
  
  
  # Calculate the standard errors and store them.
  SE <- sqrt(diag(Q))
  results$SE <- SEvec2matU(SE,k,r,p, opt)
  results$NegInvHessian <- Q
  
  
  #--------------------------------------------------------------------------------
  # GET RESIDUALS
  #--------------------------------------------------------------------------------
  
  
  epsilon <- GetResiduals(x, k, r, results$coeffs, opt)
  results$Residuals <- epsilon
  
  
  
  #--------------------------------------------------------------------------------
  # OBTAIN ROOTS OF CHARACTERISTIC POLYNOMIAL
  #--------------------------------------------------------------------------------
  
  cPolyRoots <- CharPolyRoots(results$coeffs, opt, k, r, p)
  results$cPolyRoots <- cPolyRoots
  
  
  #--------------------------------------------------------------------------------
  # PRINT OUTPUT
  #--------------------------------------------------------------------------------
  
  if (opt$print2screen) {
    
    
    if (!opt$CalcSE) {
      print(sprintf('Warning: standard errors have not been calculated!\n'))
    }
    
    # create a variable for output strings
    yesNo <- c('No','Yes')
    print(sprintf('\n-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('                      Fractionally Cointegrated VAR: Estimation Results                              '))
    print(sprintf('\n-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('Dimension of system:  %6.0f      Number of observations in sample:       %6.0f \n', p, T+opt$N))
    print(sprintf('Number of lags:       %6.0f      Number of observations for estimation:  %6.0f \n', k, T))
    print(sprintf('Restricted constant:  %6s      Initial values:                         %6.0f\n', yesNo[opt$rConstant+1], opt$N ))
    print(sprintf('Unrestricted constant:%6s      Level parameter:                        %6s\n', yesNo[opt$unrConstant+1], yesNo[opt$levelParam+1] ))
    
    
    
    
  }
  
  if (nrow(opt$R_psi) == 1) {
    
    # 1 restriction.
    dbUB <- H_psi*UB[1]
    dbLB <- H_psi*LB[1]
    dbStart <- H_psi*startVals[1]
    print(sprintf('Starting value for d:    %1.3f    Parameter space for d: (%1.3f , %1.3f) \n', dbStart[1], dbLB[1], dbUB[1]))
    print(sprintf('Starting value for b:    %1.3f    Parameter space for b: (%1.3f , %1.3f) \n', dbStart[2], dbLB[2], dbUB[2]))
  } else {
    # Unrestricted or 2 restrictions.
    print(sprintf('Starting value for d:    %1.3f    Parameter space for d: (%1.3f , %1.3f) \n', startVals[1], LB[1], UB[1]))
    print(sprintf('Starting value for b:    %1.3f    Parameter space for b: (%1.3f , %1.3f) \n', startVals[2], LB[2], UB[2]))
    print(sprintf('Imposing d >= b:      %6s\n', yesNo[opt$constrained+1] ))
    
  }
  
  print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf('Cointegrating rank:   %10.0f  AIC:            %10.3f \n', r, -2*maxLike + 2*fp))
  print(sprintf('Log-likelihood:       %10.3f  BIC:            %10.3f \n', maxLike, -2*maxLike + fp*log(T)))
  print(sprintf('log(det(Omega_hat)):  %10.3f  Free parameters:%10.0f \n', log(det(results$coeffs$OmegaHat)), fp))
  print(sprintf('-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf(    '    Fractional parameters:                                                                             \n'))
  print(sprintf(    '-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf(    '    Coefficient              \t Estimate              \t  Standard error \n'))
  print(sprintf(    '-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf(    '         d                   \t %8.3f              \t     %8.3f                \n', results$coeffs$db[1], results$SE$db[1]))
  
  if (!opt$restrictDB) {
    print(sprintf('         b                   \t %8.3f              \t     %8.3f                \n', results$coeffs$db[2], results$SE$db[2]))
  }
  
  print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
  print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
  
  
  if (r > 0) {
    
    
    if (opt$rConstant) {
      varList <- '(beta and rho):'
    }
    else {
      varList <- '(beta):        '
    }
    
    print(sprintf('    Cointegrating equations %s                                                          \n', varList))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(    '      Variable      ' ))
    
    for (j in 1:r) {
      print(sprintf(    '  CI equation %d  ', j))
    }
    
    print(sprintf('\n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    
    for (i in 1:p) {
      print(sprintf(    '        Var%d       ',i ))
      for (j in 1:r) {
        print(sprintf('    %8.3f     ', results$coeffs$betaHat[i,j] ))
      }
      print(sprintf('\n'))
    }
    
    
    if (opt$rConstant) {
      print(sprintf(    '      Constant     ' ))
      for (j in 1:r) {
        print(sprintf('    %8.3f     ', results$coeffs$rhoHat[j] ))
      }
      print(sprintf('\n'))
    }
    
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    
    if (is.null(opt$R_Alpha) & is.null(opt$R_Beta) ) {
      print(sprintf(  'Note: Identifying restriction imposed.                                                               \n'))
      print(sprintf(  '-----------------------------------------------------------------------------------------------------\n')) 
    }
    
    print(sprintf('    Adjustment matrix (alpha):                                                                         \n' ))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(    '      Variable      ' ))
    
    for (j in 1:r) {
      print(sprintf(    '  CI equation %d  ', j))
    }
    
    print(sprintf('\n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    
    for (i in 1:p) {
      print(sprintf(    '        Var %d      ',i ))
      for (j in 1:r) {
        print(sprintf('    %8.3f     ', results$coeffs$alphaHat[i,j] ))
      }
      
      print(sprintf('\n'))
      print(sprintf(    '         SE %d      ',i ))
      
      for (j in 1:r) {
        print(sprintf('   (%8.3f  )  ', results$SE$alphaHat[i,j] ))
      }
      
      print(sprintf('\n'))
    }
    
    
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(  'Note: Standard errors in parenthesis.                                                                \n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('    Long-run matrix (Pi):                                                                       \n' ))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(    '      Variable  ' ))
    
    for (j in 1:p) {
      print(sprintf(    '       Var %d   ', j))
    }
    
    print(sprintf('\n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    
    for (i in 1:p) {
      print(sprintf(    '      Var %d      ',i ))
      for (j in 1:p) {
        print(sprintf('   %8.3f    ', results$coeffs$PiHat[i,j] ))
      }
      
      print(sprintf('\n'))
    }
    
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    
    
  }
  
  
  
  
  
  # Print level parameter if present.
  if (opt$print2screen & opt$levelParam) {
    print(sprintf(  '\n-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('    Level parameter (mu):                                                                         \n' ))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    for (i in 1:p) {
      print(sprintf(    '        Var %d      ',i ))
      print(sprintf('    %8.3f     ', results$coeffs$muHat[i] ))
      print(sprintf('\n'))
      print(sprintf(    '         SE %d      ',i ))
      print(sprintf('   (%8.3f  )  ', results$SE$muHat[i] ))
      print(sprintf('\n'))
    }
    
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(  'Note: Standard errors in parenthesis (from numerical Hessian) but asymptotic distribution is unknown. \n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
  }
  
  
  
  
  # Print unrestricted constant if present.
  if (opt$print2screen & opt$unrConstant) {
    print(sprintf(  '\n-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf('    Unrestricted constant term:                                                                     \n' ))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    for (i in 1:p) {
      print(sprintf(    '        Var %d      ',i ))
      print(sprintf('    %8.3f     ', results$coeffs$xiHat[i] ))
      print(sprintf('\n'))
      print(sprintf(    '         SE %d      ',i ))
      print(sprintf('   (%8.3f  )  ', results$SE$xiHat[i] ))
      print(sprintf('\n'))
    }
    
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(  'Note: Standard errors in parenthesis (from numerical Hessian) but asymptotic distribution is unknown. \n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
  }
  
  
  
  
  # Print Gamma coefficients if required.
  if (opt$print2screen && opt$printGammas & (k > 0)) {
    
    
    for (l in 1:k) {
      
      
      GammaHatk <- results$coeffs$GammaHat[, p*(l-1)+1 : p*l]
      GammaSEk <- results$SE$GammaHat[, p*(l-1)+1 : p*l]
      
      print(sprintf('    Lag matrix %d (Gamma_%d):                                                                            \n', l, l ))
      print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
      print(sprintf(    '      Variable  ' ))
      
      for (j in 1:p) {
        print(sprintf(    '       Var %d   ', j))
      }
      
      print(sprintf('\n'))
      print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
      
      for (i in 1:p) {
        
        print(sprintf(    '      Var %d      ',i ))
        
        for (j in 1:p) {
          print(sprintf('   %8.3f    ', GammaHatk[i,j] ))
        }
        
        print(sprintf('\n'))
        print(sprintf(    '       SE %d       ',i ))
        
        for (j in 1:p) {
          print(sprintf(' (%8.3f  )  ', GammaSEk[i,j] ))
        }
        
        print(sprintf('\n'))
        
      }
      
      print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
      print(sprintf(  'Note: Standard errors in parentheses.                                                                \n'))
      print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
      
      
    }
    
    
  }
  
  
  
  
  # Print roots of characteristic polynomial if required.
  if (opt$print2screen & opt$printRoots) {
    print(sprintf(  '\n-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(  '    Roots of the characteristic polynomial                                                           \n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(  '    Number     Real part    Imaginary part       Modulus                                             \n'))
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n'))
    for (j in 1:length(cPolyRoots)) {
      print(sprintf( '      %2.0f       %8.3f       %8.3f         %8.3f                                        \n',
                     j, real(cPolyRoots[j]), imag(cPolyRoots[j]), abs(cPolyRoots[j]) ))
    }
    
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n\n'))
  }
  
  
  
  # Print notifications regarding restrictions.
  if (opt$print2screen &  (!is.null(opt$R_Alpha) | !is.null(opt$R_psi) 
                            | !is.null(opt$R_Beta) )) {
    print(sprintf(  '\n-----------------------------------------------------------------------------------------------------\n'))
    print(sprintf(  'Restrictions imposed on the following parameters:\n'))
    if(!is.null(opt$R_psi)) {
      print(sprintf('- Psi. For details see "options$R_psi"\n'))
    }
    
    if(!is.null(opt$R_Alpha)) {
      print(sprintf('- Alpha. For details see "options$R_Alpha"\n'))
    }
    
    if(!is.null(opt$R_Beta)) {
      print(sprintf('- Beta. For details see "options$R_Beta"\n'))
    }
    
    print(sprintf(  '-----------------------------------------------------------------------------------------------------\n\n'))
  }
  
  
  
  
  
  
  
  
  return(results)
}












      







################################################################################
# Define function to 
################################################################################
# 

# 
################################################################################






################################################################################
# End
################################################################################
