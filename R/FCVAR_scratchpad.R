

# Temporary staging area for Matlab code. 


function xSim <- FCVARsim(data, model, NumPeriods)
# function xSim <- FCVARsim(data, model, NumPeriods)
# Written by Michal Popiel and Morten Nielsen (This version 08.06.2015)
# 
# DESCRIPTION: This function simulates the FCVAR model as specified by
#               input "model" and starting values specified by "data."
#               Errors are drawn from a Normal distribution. 
# 
# Input <- data       (N x p matrix of data)
#         model      (a Matlab structure containing estimation results)
#         NumPeriods (number of steps for simulation)
# Output <- xSim      (NumPeriods x p matrix of simulated data values)
%_________________________________________________________________________

# --- Preliminary definitions --- %
  x <- data 
  p <- size(data,2)
  opt <- model$options
  cf  <- model$coeffs
  d <- cf$db(1)
  b <- cf$db(2)
  
  # Generate disturbance term
  err <- randn(NumPeriods,p)
  
  # --- Recursively generate simulated data values --- %
    for i <- 1:NumPeriods
  # append x with zeros to simplify calculations.
  x <- [x zeros(1,p)]
  T <- size(x,1)
  
  # Adjust by level parameter if present.
  if(opt.levelParam)
    y <- x - ones(T,1)*cf$muHat
  else
    y <- x
  end
  
  # Main term, take fractional lag.
  z <- Lbk(y,d,1)
  
  # Error correction term.
  if(~isempty(cf$alphaHat))
    z <- z + FracDiff( Lbk(y, b, 1), d-b ) * cf$PiHat'
  if(opt.rConstant)
  z <- z + FracDiff( Lbk(ones(T,1), b, 1), d-b ) * cf$rhoHat*cf$alphaHat'
  end
  end
  
  # Add unrestricted constant if present.
  if(opt.unrConstant)
    z <- z + ones(T,1)*cf$xiHat'
  end
  
  # Add lags if present.
  if(~isempty(cf$GammaHat))
  k <- size(cf$GammaHat,2)/p
  z <- z +  FracDiff(  Lbk( y , b, k)  , d) * cf$GammaHat' 
  end
  
  # Adjust by level parameter if present.
  if(opt.levelParam)
    z <- z + ones(T,1)*cf$muHat
  end
  
  # Add disturbance term
  z(T,:) <- z(T,:) + err(i,:)
  
  # Append to x matrix.
  x <- [x(1:T-1,:) z(T,:)]
  end
  # Return simulated data values, including initial values.
  xSim <- x(size(data,1)+1:end,:)
  end
  
  
  










