function [ like ] = FCVARlikeMu(y, db, mu, k, r, opt)
% function [ like ] = FCVARlikeMu(y, db, mu, k, r, opt)
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
% 
% DESCRIPTION: This function evaluates the likelihood for a given set of
% 	parameter values. It is used by the LikeGrid() function to numerically
% 	optimize over the level parameter for given values of the fractional
% 	parameters.
%
% Input = y  (matrix of variables to be included in the system)
%         db (fractional parameters d,b)
%         mu (level parameter)
%         k  (number of lags)
%         r  (number of cointegrating vectors)
%         opt (object containing the estimation options)
% Output = like (log-likelihood evaluated at specified parameter values)
%_________________________________________________________________________

    t = size(y,1);
    x = y - ones(t,1)*mu;

    % Obtain concentrated parameter estimates.
    [ estimates ] = GetParams(x, k, r, db, opt);

    % Calculate value of likelihood function.
    T = t - opt.N;
    p = size(x,2);
    like = - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates.OmegaHat));
end

