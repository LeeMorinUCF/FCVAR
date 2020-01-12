

# Temporary staging area for Matlab code. 

function [ hessian ] <- FCVARhess(x, k, r, coeffs, opt)
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



# Set dimensions of matrices.
p <- size(x, 2)

# rhoHat and betaHat are not used in the Hessian calculation.
rho <- coeffs$rhoHat
beta <- coeffs$betaHat

# Specify delta (increment for numerical derivatives).
delta <- 10^(-4)

# Convert the parameters to vector form.
phi0 <- SEmat2vecU(coeffs, k, r, p, opt)

# Calculate vector of increments.
deltaPhi <- delta*ones(size(phi0))

nPhi <- length(phi0)

# Initialize the Hessian matrix.
hessian <- zeros(nPhi)

# The loops below evaluate the likelihood at second order incremental
#   shifts in the parameters. 

for i <- 1:nPhi
for j <- 1:i

# positive shift in both parameters.
phi1 <- phi0 + deltaPhi.*( (1:nPhi) == i ) + deltaPhi.*( (1:nPhi) == j )
[ coeffsAdj ] <- SEvec2matU( phi1, k, r, p, opt )
# calculate likelihood
like1 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt )

# negative shift in first parameter, positive shift in second
# parameter. If same parameter, no shift.
phi2 <- phi0 - deltaPhi.*( (1:nPhi) == i ) + deltaPhi.*( (1:nPhi) == j )
[ coeffsAdj ] <- SEvec2matU( phi2, k, r, p, opt )
# calculate likelihood
like2 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)

# positive shift in first parameter, negative shift in second
# parameter. If same parameter, no shift.
phi3 <- phi0 + deltaPhi.*( (1:nPhi) == i ) - deltaPhi.*( (1:nPhi) == j )
[ coeffsAdj ] <- SEvec2matU( phi3, k, r, p, opt )
# calculate likelihood
like3 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)

# negative shift in both parameters.
phi4 <- phi0 - deltaPhi.*( (1:nPhi) == i ) - deltaPhi.*( (1:nPhi) == j )
[ coeffsAdj ] <- SEvec2matU( phi4, k, r, p, opt )
# calculate likelihood
like4 <- FullFCVARlike(x, k, r, coeffsAdj, beta, rho, opt)

# The numerical approximation to the second derivative.
hessian(i,j) <- ( like1 - like2 - like3 + like4 )/4/deltaPhi(i)/deltaPhi(j)

# Hessian is symmetric.
hessian(j,i) <- hessian(i,j)

end
end
end
