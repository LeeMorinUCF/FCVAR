

# Temporary staging area for Matlab code. 

function cPolyRoots <- CharPolyRoots(coeffs, opt, k, r, p)
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


b <- coeffs.db(2)

# First construct the coefficient matrix for the companion form of the VAR. 
PiStar <- eye(p)
if r > 0
PiStar <- PiStar + coeffs.alphaHat*t(coeffs.betaHat)
end

if k > 0
Gamma1 <- coeffs.GammaHat(:, 1 : p )
PiStar <- PiStar + Gamma1
for i <- 2:k

Gammai <- coeffs.GammaHat(:, (i-1)*p + 1 : i*p )
GammaiMinus1 <- coeffs.GammaHat(:, (i-2)*p + 1 : (i-1)*p )

PiStar <- [ PiStar (Gammai - GammaiMinus1) ]

end
Gammak <- coeffs.GammaHat(:, (k-1)*p + 1 : k*p )
PiStar <- [ PiStar  ( - Gammak ) ]
end

# Pad with an identity for the transition of the lagged variables.
PiStar <- [  PiStar...
            eye(p*k) zeros( p*k, p*(k>0) )  ]


# The roots are then the inverse eigenvalues of the matrix PiStar.
cPolyRoots <- 1 ./ eig(PiStar)
cPolyRoots <- sort( cPolyRoots , 'descend')

# # The inverse roots are then the eigenvalues of the matrix PiStar.
# cPolyRoots <- eig(PiStar)
# cPolyRoots <- sort( cPolyRoots , 'descend')

# Generate graph depending on the indicator plotRoots.
if (opt.plotRoots)
  # Now calculate the line for the transformed unit circle.
  # First do the negative half.
  unitCircle <- ( pi : - 0.001 : 0)
psi <- - (pi - unitCircle)/2
unitCircleX <- cos( - unitCircle)
unitCircleY <- sin( - unitCircle)
transformedUnitCircleX <- (1 - (2*cos(psi)).^b.*cos(b*psi))
transformedUnitCircleY <- (    (2*cos(psi)).^b.*sin(b*psi))
# Then do the positive half.
unitCircle <- (0 : 0.001 : pi)
psi <- (pi - unitCircle)/2
unitCircleX <- [ unitCircleX cos(unitCircle) ]
unitCircleY <- [ unitCircleY sin(unitCircle) ]
transformedUnitCircleX <- [ transformedUnitCircleX 1 ...
                           (1 - (2*cos(psi)).^b.*cos(b*psi)) ]
transformedUnitCircleY <- [ transformedUnitCircleY 0 ...
                           (    (2*cos(psi)).^b.*sin(b*psi)) ]

# Plot the unit circle and its image under the mapping
# along with the roots of the characterisitc polynomial.
figure
plot(transformedUnitCircleX, transformedUnitCircleY)
hold on
plot(unitCircleX, unitCircleY)
scatter(real(cPolyRoots), imag(cPolyRoots))
maxXYaxis <- max( [ transformedUnitCircleX unitCircleX...
                   transformedUnitCircleY unitCircleY ] )
minXYaxis <- min( [ transformedUnitCircleX unitCircleX...
                   transformedUnitCircleY unitCircleY ] )
maxXYaxis <- max( maxXYaxis, -minXYaxis )
axis( 2*[ -maxXYaxis maxXYaxis -maxXYaxis maxXYaxis ] )
hold off
axis equal
title('Roots of the characteristic polynomial with the image of the unit circle')

#     # Now transform the circle under the mapping 1/z to check if 
#     # eigenvalues are inside the circle.
#     transformedUnitCircleX1 <- transformedUnitCircleX./ ... 
#                     (transformedUnitCircleX.^2+transformedUnitCircleY.^2)
#     transformedUnitCircleY1 <- -transformedUnitCircleY./ ...
#                     (transformedUnitCircleX.^2+transformedUnitCircleY.^2)

#     # Plot the unit circle and its image under the mapping
#     # along with the inverse roots of the characterisitc polynomial.
#     plot(transformedUnitCircleX1, transformedUnitCircleY1)
#     hold on
#     plot(unitCircleX, unitCircleY)
#     scatter(real(cPolyRoots), imag(cPolyRoots))
#     maxXYaxis <- max( [ transformedUnitCircleX1 unitCircleX...
#                        transformedUnitCircleY1 unitCircleY ] )
#     minXYaxis <- min( [ transformedUnitCircleX1 unitCircleX...
#                        transformedUnitCircleY1 unitCircleY ] )
#     maxXYaxis <- max( maxXYaxis, -minXYaxis )
#     axis( 2*[ -maxXYaxis maxXYaxis -maxXYaxis maxXYaxis ] )
#     hold off
#     title('Inverse roots of the characteristic polynomial with the image of the unit circle')

end

end


