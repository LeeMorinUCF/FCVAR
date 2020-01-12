

# Temporary staging area for Matlab code. 


function [ param ] <- SEmat2vecU( coeffs, k, r, p , opt)
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
%_________________________________________________________________________

# If restriction d=b is imposed, only adjust d.
if opt$restrictDB
param <- coeffs$db(1)
else
  param <- coeffs$db
end

# Level parameter MuHat.
if opt$levelParam
param <- [ param reshape(coeffs$muHat,1,p) ]
end

# Unrestricted constant.
if opt$unrConstant
param <- [ param reshape(coeffs$xiHat,1,p)]
end

# alphaHat
if r > 0
param <- [  param reshape( coeffs$alphaHat, 1, p*r ) ]
end

# GammaHat
if k > 0
param <- [  param reshape( coeffs$GammaHat, 1, p*p*k ) ]
end

end




function [ coeffs ] <- SEvec2matU( param, k, r, p, opt )
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
%_________________________________________________________________________

if opt$restrictDB
coeffs$db <- [ param(1) param(1) ]  # store d,b
param <- param(2:end)				# drop d,b from param
else
  coeffs$db <- param(1:2)
param <- param(3:end)
end

if opt$levelParam
coeffs$muHat <- param(1:p)
param <- param(p+1:end)
end

if opt$unrConstant
coeffs$xiHat <- matrix(param(1:p), nrow = p, ncol = 1)
param <- param(p+1:end)
end

if r > 0
coeffs$alphaHat <- reshape( param(1:p*r), p, r)
param <- param(p*r+1:end)
else
coeffs$alphaHat <- NULL
end

if k > 0
coeffs$GammaHat <- reshape( param(1 : end), p, p*k)
else
coeffs$GammaHat <- NULL
end

end


