function [ hessian ] = FCVARhess(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
% 
% function [ hessian ] = FCVARhess(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% FCVARhess(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions) is the Hessian matrix
%       for the fractionally cointegrated VAR model evaluated at {db, alpha, beta, rho, Gamma}.
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
% 
% output = matrix hessian is the Hessian matrix for the fractionally
%       cointegrated VAR model evaluated at {db, alpha, beta, rho, Gamma}.
% 
% 
%
% Calls the functions FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions)
%       to get the likelihood and also calls SEmat2vec( db, alpha, Gamma, k, r, p )
%       and SEvec2mat( param, k, r, p ) to convert parameter vector to matrices and
%       vice versa.
% 
%______________________________________________________


% Set dimensions of matrices.
p = size(x, 2);

% Specify delta.
delta = 10^(-4);

% Convert the parameters to vector form.
phi0 = SEmat2vec( db, alpha, Gamma, k, r, p );

% Calculate vector of increments.
deltaPhi = delta*ones(size(phi0));

nPhi = length(phi0);
for i = 1:nPhi
    
    
    for j = 1:i
        
        phi1 = phi0 + deltaPhi.*( (1:nPhi) == i ) + deltaPhi.*( (1:nPhi) == j );
        [ db, alpha, Gamma ] = SEvec2mat( phi1, k, r, p );
        like1 = FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions);
        
        phi2 = phi0 - deltaPhi.*( (1:nPhi) == i ) + deltaPhi.*( (1:nPhi) == j );
        [ db, alpha, Gamma ] = SEvec2mat( phi2, k, r, p );
        like2 = FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions);
        
        phi3 = phi0 + deltaPhi.*( (1:nPhi) == i ) - deltaPhi.*( (1:nPhi) == j );
        [ db, alpha, Gamma ] = SEvec2mat( phi3, k, r, p );
        like3 = FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions);
        
        phi4 = phi0 - deltaPhi.*( (1:nPhi) == i ) - deltaPhi.*( (1:nPhi) == j );
        [ db, alpha, Gamma ] = SEvec2mat( phi4, k, r, p );
        like4 = FullFCVARlike(x, k, r, db, alpha, beta, rho, Gamma, EstnOptions);
        
        hessian(i,j) = ( like1 - like2 - like3 + like4 )/4/deltaPhi(i)/deltaPhi(j);
        hessian(j,i) = hessian(i,j);
    end
end




% end


