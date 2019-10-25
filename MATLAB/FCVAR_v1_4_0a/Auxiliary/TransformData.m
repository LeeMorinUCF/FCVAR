function [ Z0, Z1, Z2, Z3 ] = TransformData(x, k, db, opt)
% function [ Z0, Z1, Z2, Z3 ] = TransformData(x, k, db, opt)
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
% Based on Lee Morin & Morten Nielsen (May 24, 2013)
%
% DESCRIPTION: Returns the transformed data required for regression and
% 	reduced rank regression.
% 
% Input = x   (matrix of variables to be included in the system)
%         k   (number of lags)
%         db  (fractional differencing parameters d and b)
%         opt (object containing the estimation options)
% Output = Z0, Z1, Z2, and Z3 of transformed data.
%
% Calls the function FracDiff(x, d) and Lbk(x, b, k).
% _____________________________________________________
 
	% Number of initial values and sample size.
    N = opt.N;
    T = size(x,1) - N;

    % Extract parameters from input.
    d = db(1);
    b = db(2);

    % Transform data as required.
    Z0 = FracDiff(x, d);
    
    Z1 = x;
    % Add a column with ones if model includes a restricted constant term.
    if(opt.rConstant)
        Z1 = [ x ones(N+T,1) ];
    end
    Z1 = FracDiff(  Lbk( Z1 , b, 1)  ,  d - b );
    
    Z2 = FracDiff(  Lbk( x , b, k)  , d);

    % Z3 contains the unrestricted deterministics
    Z3 = [];
    % % Add a column with ones if model includes a unrestricted constant term.
    if(opt.unrConstant)
        Z3 = ones(T,1);
    end

    % Trim off initial values.
    Z0 = Z0( N+1:end ,:);
    Z1 = Z1( N+1:end ,:);
    if(k>0)
        Z2 = Z2( N+1:end ,:);
    end
    
end
