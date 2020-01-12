

# Temporary staging area for Matlab code. 



function [ fp ] <- FreeParams(k, r, p, opt, rankJ)
# function [ fp ] <- FreeParams(k, r, p, opt, rankJ)
# Written by Michal Popiel and Morten Nielsen (This version 05.22.2015)
# 
# DESCRIPTION: This function counts the number of free parameters based on
# 	the number of coefficients to estimate minus the total number of
# 	restrictions. When both alpha and beta are restricted, the rank condition
# 	is used to count the free parameters in those two variables.
%
# Input <- x (matrix of variables to be included in the system)
#         k (number of lags)
#         r (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- fp (number of free parameters)
%_________________________________________________________________________


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
    rDB <- size(opt$R_psi,1)
  
  
  if(isempty(opt$R_Beta) && isempty(opt$R_Alpha))
    # If Alpha or Beta are unrestricted, only an identification restriction
  #   is imposed
  rB <- r*r # identification restrictions (eye(r))
  numRest <- rDB + rB
  # ------ Free parameters is the difference between the two ---- %
    fp <- numParams - numRest
  else
    # If there are restrictions imposed on beta or alpha, we can use the rank
  #   condition calculated in the main estimation function to give us the
  #   number of free parameters in alpha, beta, and the restricted constant
  #   if it is being estimated:
    fpAlphaBetaRhoR <- rankJ
  fp <- fpAlphaBetaRhoR + (fDB + fpG + fpM + fpUrh) - rDB    
  end
  
  end
  
