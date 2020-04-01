################################################################################
# 
# Library of Functions for FCVAR Estimation
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
# DESCRIPTION: A library of functions called by the FCVAR estimation function. 
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
# Define function to calculate cointegrating parameters from the data. 
################################################################################
# 
# function [ estimates ] <- GetParams(x, k, r, db, opt)
# Written by Michal Popiel and Morten Nielsen (This version 04.13.2016)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function uses FWL and reduced rank regression to obtain
# the estimates of Alpha, Beta, Rho, Pi, Gamma, and Omega
# 
# Input <- x   (matrix of variables to be included in the system)
#         k   (number of lags)
#         r   (number of cointegrating vectors)
#         db  (value of d and b)
#         opt (object containing the estimation options)
# Output <- estimates (Matlab structure containing the following)
#          - estimates.db (taken directly from the input)
#          - estimates.alphaHat
#          - estimates.betaHat
#          - estimates.rhoHat
#          - estimates.piHat
#          - estimates.OmegaHat
#          - estimates.GammaHat ( p x kp matrix [GammaHat1,...,GammaHatk])
#
################################################################################


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




################################################################################
# Define function to calculate the likelihood function
# on a grid of candidate parameter values.
################################################################################
# 
# function [ params ] <- LikeGrid(x,k,r,opt)
# Written by Michal Popiel and Morten Nielsen (This version 04.12.2016)
# 
# DESCRIPTION: This function evaluates the likelihood over a grid of values
# 	for (d,b) (or phi). It can be used when parameter estimates are sensitive to
# 	starting values to give an approximation of the global max which can
# 	then be used as the starting value in the numerical optimization in
# 	FCVARestn().
# 
# Input <- x   (matrix of variables to be included in the system)
#         k   (number of lags)
#         r   (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- params (row vector of d,b, and mu (if level parameter is selected)
#					corresponding to a maximum over the grid of (d,b), or phi)
# 
# Note:	If opt$LocalMax == 0, LikeGrid returns the parameter values
#       corresponding to the global maximum of the likelihood on the grid.
#       If opt$LocalMax == 1, LikeGrid returns the parameter values for the
#       local maximum corresponding to the highest value of b. This
#       alleviates the identification problem mentioned in Johansen and
#       Nielsen (2010, section 2.3).
# 
################################################################################


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










################################################################################
# Define function to calculate the likelihood for the unconstrained model. 
################################################################################
# 
# function [ like ] <- FCVARlikeMu(y, db, mu, k, r, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# 
# DESCRIPTION: This function evaluates the likelihood for a given set of
# 	parameter values. It is used by the LikeGrid() function to numerically
# 	optimize over the level parameter for given values of the fractional
# 	parameters.
# 
# Input <- y  (matrix of variables to be included in the system)
#         db (fractional parameters d,b)
#         mu (level parameter)
#         k  (number of lags)
#         r  (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- like (log-likelihood evaluated at specified parameter values)
# 
################################################################################



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




################################################################################
# Define function to calculate the likelihood for the constrained model. 
################################################################################
# 
# function [ like ] <- FCVARlike(x, params, k, r, opt)
# Written by Michal Popiel and Morten Nielsen (This version 11.10.2014)
# 
# DESCRIPTION: This function adjusts the variables with the level parameter,
# 	if present, and returns the log-likelihood given d,b.
# 
# Input <- x   (matrix of variables to be included in the system)
#         params (a vector of parameters d,b, and mu (if option selected))
#         k   (number of lags)
#         r   (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- like (concentrated log-likelihood evaluated at given parameters)
# 
################################################################################



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




################################################################################
# Define function to calculate the likelihood, 
# evaluated at the parameters provided as inputs. 
################################################################################
# 
# function [ like ] <- FullFCVARlike(x, k, r, coeffs, beta, rho, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function returns the value of the log-likelihood
# 	evaluated at the parameters provided as inputs.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         coeffs (Matlab structure of coefficients)
#         beta (value of beta)
#         rho  (value of rho)
#         opt  (object containing the estimation options)
# Output <- like (value of the log likelihood)
# 
################################################################################

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



################################################################################
# Define function to transform the raw data by fractional differencing.
################################################################################
# 
# function [ Z0, Z1, Z2, Z3 ] <- TransformData(x, k, db, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (May 24, 2013)
# 
# DESCRIPTION: Returns the transformed data required for regression and
# 	reduced rank regression.
# 
# Input <- x   (matrix of variables to be included in the system)
#         k   (number of lags)
#         db  (fractional differencing parameters d and b)
#         opt (object containing the estimation options)
# Output <- Z0, Z1, Z2, and Z3 of transformed data.
# 
# Calls the function FracDiff(x, d) and Lbk(x, b, k).
# 
################################################################################

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


################################################################################
# Define function to calculate restricted estimates with a switching algorithm.
################################################################################
# 
# function [betaStar, alphaHat, OmegaHat] 
#                   <- RstrctOptm_Switch(beta0, S00, S01, S11, T, p, opt)
# Written by Michal Popiel and Morten Nielsen (This version 03.29.2016)
# 
# DESCRIPTION: This function is imposes the switching algorithm of Boswijk
#   and Doornik (2004, page 455) to optimize over free parameters psi 
#   and phi directly, combined with the line search proposed by 
#	Doornik (2016, working paper). We translate between  (psi, phi) and 
#	(alpha, beta) using the relation of R_Alpha*vec(alpha) <- 0 and 
#	A*psi <- vec(alpha'), and R_Beta*vec(beta) <- r_beta and 
#	H*phi+h <- vec(beta). Note the transposes.
#
# Input <- beta0 (unrestricted estimate of beta)
#         S00, S01, S11 (product moments)
#         T (number of observations)
#         p (number of variables)
#         opt (object containing the estimation options)
# Output <- betaStar (estimate of betaStar)
#          alphaHat (estimate of alpha)
#          OmegaHat (estimate of Omega)
# 
################################################################################


RstrctOptm_Switch <- function(beta0, S00, S01, S11, T, p, opt) {
  
  
  r  <- ncol(beta0)
  p1 <- p + opt$rConstant
  
  # Restrictions on beta.
  if(is.null(opt$R_Beta)) {
    H <- diag(p1*r)
    h <- matrix(0, nrow = p1*r, ncol = 1)
  } else {
    H <- null(opt$R_Beta)
    h <- t(opt$R_Beta) %*% 
      solve(opt$R_Beta %*% t(opt$R_Beta)) %*% opt$r_Beta 
  }
  
  
  # Restrictions on alpha.
  #   We use the commutation matrix K_pr to transform vec(A) into vec(A'),
  #   see Magnus & Neudecker (1988, p. 47, eqn (1)).
  Ip  <- diag(p)
  Kpr <- matrix(kron(Ip, diag(r)), nrow = p*r, ncol = p*r)
  if(is.null(opt$R_Alpha)) {
    A <- Kpr %*% diag(p*r)  
  } else {
    A <- null(opt$R_Alpha %*% solve(Kpr))
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
  # print('kron(alphaHat, diag(p1)) = ')
  # print(kron(alphaHat, diag(p1)))
  # print('h = ')
  # print(h)
  # print('kron(alphaHat, diag(p1)) %*% h = ')
  # print(kron(alphaHat, diag(p1)) %*% h)
  # 
  
  
  
  # Get candidate values for entering the switching algorithm.
  vecPhi1 <- solve(t(H) %*% 
                     kron(t(alphaHat) %*% solve(OmegaHat) %*% alphaHat, S11) %*% 
                     H) %*% 
    t(H) %*% (kron(t(alphaHat) %*% solve(OmegaHat), S11)) %*% 
    (vecPiLS - kron(alphaHat, diag(p1)) %*% h)
  
  # Translate vecPhi to betaStar.
  vecB <- H %*% vecPhi1 + h
  betaStar <- matrix(vecB, nrow = p1, ncol = r)    
  
  # Candidate value of vecPsi.
  vecPsi1 <- solve(t(A) %*% 
                     kron(solve(OmegaHat), t(betaStar) %*% S11 %*% betaStar) %*% 
                     A) %*%
    t(A) %*% (kron(solve(OmegaHat), t(betaStar) %*% S11)) %*% vecPiLS
  
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
    vecPsi1 <- solve(t(A) %*% kron(solve(OmegaHat), t(betaStar) %*% 
                                     S11 %*% betaStar) %*% A) %*%
      t(A) %*% (kron(solve(OmegaHat), t(betaStar) %*% S11)) %*% vecPiLS
    
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
                       kron(t(alphaHat) %*% inv(OmegaHat) %*% alphaHat, S11) %*% 
                       H) %*% 
      t(H) %*% (kron(t(alphaHat) %*% solve(OmegaHat), S11)) %*%
      (vecPiLS - kron(alphaHat, diag(p1)) %*% h)
    
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



################################################################################
# Define function to calculate residuals from given parameter values.
################################################################################
# 
# function [ epsilon ] <- GetResiduals(x, k, r, coeffs, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function calculates the model residuals.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         coeffs (Matlab structure of coefficients)
#         opt (object containing the estimation options)
# Output <- epsilon (matrix of residuals from model estimation evaluated at
#                     the parameter estimates specified in coeffs) 
# 
################################################################################

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





################################################################################
# Define function to calculate lag polynomial in the fractional lag operator.
################################################################################
# 
# function [ Lbkx ] <- Lbk(x, b, k)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (May 24, 2013)
# 
# DESCRIPTION: Lbk(x, b, k) is a lag polynomial in the fractional lag operator.
# 
# Input <- x (vector or matrix of data)
#         b (scalar value at which to calculate the fractional lag)
#         k (number of lags)
# 
# Output <- matrix  [ Lb^1 x, Lb^2 x, ..., Lb^k x] where Lb <- 1 - (1-L)^b. The output 
#		matrix has the same number of rows as x but k times as many columns.
# 
# Calls the function FracDiff(x, d) 
# 
################################################################################


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



################################################################################
# Define function to fractionally difference a series. 
################################################################################
# 
# function [dx] = FracDiff(x,d)
# Lee Morin & Morten Nielsen
# 
# 
# FracDiff(x,d) is a fractional differencing procedure based on the
# 	simple filter in the time domain.
#   It is used mainly to test other calculations. 
# 
# input = x (vector or matrix of data)
#         d (scalar value at which to calculate the fractional difference)
# 
# output = vector or matrix (1-L)^d x of same dimensions as x.
# 
################################################################################

FracDiff_filter <- function(x, d) {
  
  if(is.null(x)) {
    dx <- NULL
  } else {
    
    
    T <- nrow(x)
    p <- ncol(x)
    
    # print('T = ')
    # print(T)
    # print('p = ')
    # print(p)
    
    # print(summary(x))
    
    k <- matrix(seq(1, T-1), nrow = T-1, ncol = 1)
    b <- (k-d-1)/k
    b <- c(1, cumprod(b))
    
    # xlead <- data.frame(matrix(0,
    #                            nrow = T,
    #                            ncol = ncol(x)))
    xlead <- matrix(0,
                    nrow = T,
                    ncol = ncol(x))
    colnames(xlead) <- colnames(x)
    
    # print('xlead = ')
    # print(xlead)
    # print('x = ')
    # print(x)
    # print('rbind(xlead, x) = ')
    # print(rbind(xlead, x))
    
    
    
    dx <- filter(rbind(xlead, x),
                 filter = b,
                 sides = 1)
    # Trim off the leading rows.
    dx <- dx[seq(T+1, 2*T), ]
    
  }
  
  return(dx)
}


################################################################################
# Define function to fractionally difference a series. 
################################################################################
# 
# function [dx] = FracDiff(x,d)
# Andreas Noack Jensen & Morten Nielsen
# May 24, 2013
# 
# FracDiff(x,d) is a fractional differencing procedure based on the
# 	fast fractional difference algorithm of Jensen & Nielsen (2014, JTSA).
# 
# input = x (vector or matrix of data)
#         d (scalar value at which to calculate the fractional difference)
# 
# output = vector or matrix (1-L)^d x of same dimensions as x.
# 
################################################################################


FracDiff <- function(x, d) {
  
  
  # print(summary(x))
  
  if(is.null(x)) {
    dx <- NULL
  } else {
    
    
    # T <- nrow(x)
    p <- ncol(x)
    
    
    
    iT <- nrow(x)
    np2 <- nextn(2*iT - 1, 2)
    
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
      
      
      dxi <- fft(fft(c(b, rep(0, np2 - iT))) * fft(c(x[, i], rep(0, np2 - iT))), inverse = T) / np2
      
      dx[, i] <- Re(dxi[1:iT])
      
    }
    
    
    
  }
  
  return(dx)
}



################################################################################
# Define function to count the number of free parameters
################################################################################
# 
# function [ fp ] <- FreeParams(k, r, p, opt, rankJ)
# Written by Michal Popiel and Morten Nielsen (This version 05.22.2015)
# 
# DESCRIPTION: This function counts the number of free parameters based on
# 	the number of coefficients to estimate minus the total number of
# 	restrictions. When both alpha and beta are restricted, the rank condition
# 	is used to count the free parameters in those two variables.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- fp (number of free parameters)
# 
################################################################################


FreeParams <- function(k, r, p, opt, rankJ) {
  
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
  
  
  if(isempty(opt$R_Beta) & isempty(opt$R_Alpha)) {
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


################################################################################
# Define function to calculate the Hessian matrix. 
################################################################################
# 
# function [ hessian ] <- FCVARhess(x, k, r, coeffs, opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# 
# DESCRIPTION: This function calculates the Hessian matrix of the
# 	log-likelihood numerically.
# 
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         coeffs (coefficient estimates around which estimation takes place)
#         opt (object containing the estimation options)
# Output <- hessian (matrix of second derivatives)
# 
################################################################################



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
      like1 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt )
      
      # negative shift in first parameter, positive shift in second
      # parameter. If same parameter, no shift.
      phi2 <- phi0 - deltaPhi*( (1:nPhi) == i ) + deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi2, k, r, p, opt )
      # calculate likelihood
      like2 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)
      
      # positive shift in first parameter, negative shift in second
      # parameter. If same parameter, no shift.
      phi3 <- phi0 + deltaPhi*( (1:nPhi) == i ) - deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi3, k, r, p, opt )
      # calculate likelihood
      like3 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)
      
      # negative shift in both parameters.
      phi4 <- phi0 - deltaPhi*( (1:nPhi) == i ) - deltaPhi*( (1:nPhi) == j )
      coeffsAdj <- SEvec2matU( phi4, k, r, p, opt )
      # calculate likelihood
      like4 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)
      
      # The numerical approximation to the second derivative.
      hessian[i,j] <- ( like1 - like2 - like3 + like4 )/4/deltaPhi[i]/deltaPhi[j]
      
      # Hessian is symmetric.
      hessian[j,i] <- hessian[i,j]
      
      
    }
    
  }
  
  return(hessian)
}


################################################################################
# Define function to collect parameters into a vector
################################################################################
# 
# function [ param ] <- SEmat2vecU( coeffs, k, r, p , opt)
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function transforms the model parameters in matrix
# 	form into a vector.
# 
# Input <- coeffs (Matlab structure of coefficients in their usual matrix form)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         p (number of variables in the system)
#         opt (object containing the estimation options)
# Output <- param (vector of parameters)
# 
################################################################################


SEmat2vecU <- function( coeffs, k, r, p , opt) {
  
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



################################################################################
# Define function to extract parameters from a vector. 
################################################################################
# 
# function [ coeffs ] <- SEvec2matU( param, k, r, p, opt )
# Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
# Based on Lee Morin & Morten Nielsen (August 22, 2011)
# 
# DESCRIPTION: This function transforms the vectorized model parameters
# 	into matrices.
# 
# Input <- param (vector of parameters)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         p (number of variables in the system)
#         opt (object containing the estimation options)
# Output <- coeffs (Matlab structure of coefficients in their usual matrix form)
# 
################################################################################

SEvec2matU <- function( param, k, r, p, opt ) {
  
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




################################################################################
# Define function to calculate the roots of the characteristic polynomial
################################################################################
# 
# function cPolyRoots <- CharPolyRoots(coeffs, opt, k, r, p)
# Written by Michal Popiel and Morten Nielsen (This version 12.07.2015)
# Based on Lee Morin & Morten Nielsen (May 31, 2013)
# 
# DESCRIPTION: CharPolyRoots calculates the roots of the 
#     characteristic polynomial and plots them with the unit circle 
#     transformed for the fractional model, see Johansen (2008).
# 
# input <- coeffs (Matlab structure of coefficients
#         opt (object containing the estimation options)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         p (number of variables in the system)
# 
# output <- complex vector cPolyRoots with the roots of the characteristic polynomial.
# 
# No dependencies.
# 
# Note: The roots are calculated from the companion form of the VAR, 
#       where the roots are given as the inverse eigenvalues of the 
#       coefficient matrix.
# 
################################################################################

CharPolyRoots <- function(coeffs, opt, k, r, p) {
  
  
  b <- coeffs$db[2]
  
  # First construct the coefficient matrix for the companion form of the VAR. 
  PiStar <- diag(p)
  if (r > 0) {
    PiStar <- PiStar + coeffs$alphaHat %*% t(coeffs$betaHat)
  }
  
  # print('coeffs$GammaHat = ')
  # print(coeffs$GammaHat)
  # print('i = ')
  # print(i)
  # print('p = ')
  # print(p)
  # print('2:k = ')
  # print(2:k)
  
  
  if (k > 0) {
    Gamma1 <- coeffs$GammaHat[ , 1 : p]
    PiStar <- PiStar + Gamma1
    if (k > 1) {
      for (i in 2:k) {
        
        Gammai <- coeffs$GammaHat[ , seq(((i-1)*p + 1), i*p)]
        GammaiMinus1 <- coeffs$GammaHat[ , seq(((i-2)*p + 1), (i-1)*p)]
        
        PiStar <- cbind(PiStar, (Gammai - GammaiMinus1))
        
      }
    }
    
    Gammak <- coeffs$GammaHat[ , seq(((k-1)*p + 1), k*p)]
    PiStar <- cbind(PiStar, ( - Gammak ))
  }
  
  # print('PiStar = ')
  # print(PiStar)
  # print('p = ')
  # print(p)
  # print('k = ')
  # print(k)
  
  # Pad with an identity for the transition of the lagged variables.
  if (k > 0) {
    PiStar <- rbind(PiStar, 
                    cbind(diag(p*k), 
                          matrix(0, nrow = p*k, ncol = p )))
  }
  
   
  # The roots are then the inverse eigenvalues of the matrix PiStar.
  cPolyRoots <- 1 / eigen(PiStar)$values
  cPolyRoots <- cPolyRoots[order(-Mod(cPolyRoots))]
  
  # Generate graph depending on the indicator plotRoots.
  if (opt$plotRoots) {
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
    
    # Determine axes based on largest roots. 
    maxXYaxis <- max( c(transformedUnitCircleX, unitCircleX,
                        transformedUnitCircleY, unitCircleY) )
    minXYaxis <- min( c(transformedUnitCircleX, unitCircleX,
                        transformedUnitCircleY, unitCircleY) )
    maxXYaxis <- max( maxXYaxis, -minXYaxis )
    
    plot(transformedUnitCircleX, 
         transformedUnitCircleY, 
         main = c('Roots of the characteristic polynomial', 
                  'with the image of the unit circle'), 
         xlab = 'Real Part of Root', 
         ylab = 'Imaginary Part of Root', 
         xlim = 2*c(-maxXYaxis, maxXYaxis),
         ylim = 2*c(-maxXYaxis, maxXYaxis),
         type = 'l', 
         lwd = 3, 
         col = 'red')
    lines(unitCircleX, unitCircleY, lwd = 3, col = 'black')
    points(Re(cPolyRoots), Im(cPolyRoots), 
           pch = 16, col = 'blue')
    
    
    
  }
  
  
  return(cPolyRoots)
}






# Functions to make: 


#--------------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------------



################################################################################
# Define function to 
################################################################################
# 

# 
################################################################################




################################################################################
# End
################################################################################

