function [ like ] = FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
% 
% function like = FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions) is the unconcentrated likelihood 
%       function for the fractionally cointegrated VAR model.
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
% output = scalar like is the value of the log likelihood function evaluated at 
%       the specified parameters. 
% 
% Calls the functions GetResiduals(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
%       to get the residuals.
% 
%______________________________________________________

% Determine desired set of options.
N = EstnOptions{1};

% Obtain residuals.
[ epsilon ] = GetResiduals(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions);


% Calculate value of likelihood function.
T = size(x,1) - N;
p = size(x,2);
OmegaHat = epsilon'*epsilon/T;
like = - T*p/2*( log(2*pi) + 1)  - T/2*log(det(OmegaHat));



% end


