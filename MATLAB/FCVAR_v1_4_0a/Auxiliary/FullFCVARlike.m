function [ like ] = FullFCVARlike(x, k, r, coeffs, beta, rho, opt)
% function [ like ] = FullFCVARlike(x, k, r, coeffs, beta, rho, opt)
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
% Based on Lee Morin & Morten Nielsen (August 22, 2011)
% 
% DESCRIPTION: This function returns the value of the log-likelihood
% 	evaluated at the parameters provided as inputs.
%
% Input = x (matrix of variables to be included in the system)
%         k (number of lags)
%         r (number of cointegrating vectors)
%         coeffs (Matlab structure of coefficients)
%         beta (value of beta)
%         rho  (value of rho)
%         opt  (object containing the estimation options)
% Output = like (value of the log likelihood)
%_________________________________________________________________________

% Add betaHat and rhoHat to the coefficients to get residuals because they
%   are not used in the Hessian calculation and are missing from the
%   structure coeffs
coeffs.betaHat = beta;
coeffs.rhoHat = rho;

% Obtain residuals.
[ epsilon ] = GetResiduals(x, k, r, coeffs, opt);


% Calculate value of likelihood function.
T = size(x,1) - opt.N;
p = size(x,2);
OmegaHat = epsilon'*epsilon/T;
like = - T*p/2*( log(2*pi) + 1)  - T/2*log(det(OmegaHat));

end
