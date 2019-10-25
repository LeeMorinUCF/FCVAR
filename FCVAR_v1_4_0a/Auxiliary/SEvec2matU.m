function [ coeffs ] = SEvec2matU( param, k, r, p, opt )
% function [ coeffs ] = SEvec2matU( param, k, r, p, opt )
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
% Based on Lee Morin & Morten Nielsen (August 22, 2011)
% 
% DESCRIPTION: This function transforms the vectorized model parameters
% 	into matrices.
%
% Input = param (vector of parameters)
%         k (number of lags)
%         r (number of cointegrating vectors)
%         p (number of variables in the system)
%         opt (object containing the estimation options)
% Output = coeffs (Matlab structure of coefficients in their usual matrix form)
%_________________________________________________________________________

    if opt.restrictDB
        coeffs.db = [ param(1) param(1) ];  % store d,b
        param = param(2:end);				% drop d,b from param
    else
        coeffs.db = param(1:2);
        param = param(3:end);
    end

    if opt.levelParam
        coeffs.muHat = param(1:p);
        param = param(p+1:end);
    end

    if opt.unrConstant
        coeffs.xiHat = param(1:p)';
        param = param(p+1:end);
    end
    
    if r > 0
        coeffs.alphaHat = reshape( param(1:p*r), p, r);
        param = param(p*r+1:end);
    else
        coeffs.alphaHat = [];
    end

    if k > 0
        coeffs.GammaHat = reshape( param(1 : end), p, p*k);
    else
        coeffs.GammaHat = [];
    end

end

