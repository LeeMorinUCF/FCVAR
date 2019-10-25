function [ like ] = FCVARlike(x, k, r, db, EstnOptions)
% 
% function like = FCVARlike(x, k, r, db, options)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% FCVARlike(x, k, r, db, EstnOptions) is the likelihood function for the fractional AR model.
% 
% input = vector or matrix of data x.
%       scalar k denotes the number of lags.
%       scalar r denotes the cointegrating rank.
%       vector or scalar db of fractional differencing parameters.
%           The size of db indicates whether or not restriction d = b is imposed
%           based on the number of free parameters.
%       cell array EstnOptions.
% 
%   EstnOptions is an array of options for estimating the remaining parameters.
%   EstnOptions = { N, M, deterministics, R_Gamma, r_Gamma, [], [], [],
%                   nCols, A, H, a, b, h, [], [], [], [],
%                   [], [], [], print2screen, printGammas,
%                   restrictFminOptions }
% The elements of EstnOptions are described in detail in the comments in
%       functions DefaultEstnOptions(db0, k, r, p),
%       and RestrictEstnOptions(db0, k, r, p, defaultFCVARoptions)
%       and AdjustEstnOptions(db0, k, r, p, initialFCVARoptions).
% 
% output = scalar like is the value of the log likelihood function evaluated at db, 
%       with all other parameters concentrated out. 
% 
% 
%
% Calls the functions GetParams(x, k, r, db, EstnOptions) to get the OmegaHat matrix.
% 
%______________________________________________________


% Determine desired set of options.
N = EstnOptions{1};

% Obtain concentrated parameter estimates.
[ estimates ] = GetParams(x, k, r, db, EstnOptions);
% estimates = { alphaHat, betaHat, rhoHat, PiHat, OmegaHat, GammaHat } 
OmegaHat = cell2mat(estimates(5));

% Calculate value of likelihood function.
T = size(x,1) - N;
p = size(x,2);
like = - T*p/2*( log(2*pi) + 1)  - T/2*log(det(OmegaHat));



% end


