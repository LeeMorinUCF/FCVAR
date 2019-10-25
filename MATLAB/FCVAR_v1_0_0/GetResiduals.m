function [ epsilon ] = GetResiduals(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
% 
% function [ epsilon ] = GetResiduals(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% GetResiduals(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions) obtains the residuals 
%   for the fractional AR model.
% 
% input = vector or matrix of data x.
%       scalar k denotes the number of lags.
%       scalar r denotes the cointegrating rank.
%       vector or scalar db of fractional differencing parameters.
%           The size of db indicates whether or not restriction d = b is imposed
%           based on the number of free parameters.
%       matrices alpha, beta, rho, Gamma are parameter values for the FCVAR model.
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
% output = matrix of residuals from the fractionally cointegrated VAR model 
%       evaluated at the parameter values in the inputs. 
% 
% 
%
% Calls the function TransformData(x, k, db, EstnOptions) to transform the data.
% 
%______________________________________________________


% Transform data.
[ Z0, Z1, Z2 ] = TransformData(x, k, db, EstnOptions);


% Calculate residuals.
epsilon = Z0;

if r > 0
    epsilon = epsilon - Z1*[ beta; rho ]*alpha';
end

if k > 0
    epsilon = epsilon - Z2*Gamma';
end

  

% end


