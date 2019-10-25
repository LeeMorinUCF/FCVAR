function [ Z0, Z1, Z2 ] = TransformData(x, k, db, EstnOptions)

% 
% function [ Z0, Z1, Z2 ] = TransformData(x, k, db, EstnOptions)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% 
% TransformData(x, k, db, EstnOptions) returns the transformed data required for 
%       regression and reduced rank regression.
% 
% 
% input = vector or matrix of data x.
%       scalar k denotes the number of lags.
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
% output = matrices of transformed data.
%
% Calls the functions FracDiff(x, d, M) and Lbk(x, b, k, M).
% _____________________________________________________


% Determine desired set of options.
N = EstnOptions{1};
M = EstnOptions{2};
deterministics = EstnOptions(3);

% Extract parameters from input.

d = db(1);
% Restriction d = b imposed if d is a scalar.
if length(db) == 1
    b = d;
else
    b = db(2);
end
T = size(x, 1) - N;
p = size(x, 2);


% Transform data as required.
Z0 = FracDiff(x, d, M);
Z2 = FracDiff(  Lbk( x , b, k, M)  , d, M);
Z1 = x;
% Add a column of ones if model includes a restricted constant term.
if strcmp(deterministics, 'restricted constant')
    % For restriced constant, this includes a column of ones.
    Z1 = [ x [ zeros(N,1); ones(T,1) ] ];
end
Z1 = FracDiff(  Lbk( Z1 , b, 1, M)  ,  d - b, M );


% Trim off initial values.
Z0 = Z0( N+1:end ,:);
Z2 = Z2( N+1:end ,:);
Z1 = Z1( N+1:end ,:);






% end


