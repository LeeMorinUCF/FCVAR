function [ param ] = SEmat2vecU( coeffs, k, r, p , opt)
% function [ param ] = SEmat2vecU( coeffs, k, r, p , opt)
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
% Based on Lee Morin & Morten Nielsen (August 22, 2011)
% 
% DESCRIPTION: This function transforms the model parameters in matrix
% 	form into a vector.
%
% Input = coeffs (Matlab structure of coefficients in their usual matrix form)
%         k (number of lags)
%         r (number of cointegrating vectors)
%         p (number of variables in the system)
%         opt (object containing the estimation options)
% Output = param (vector of parameters)
%_________________________________________________________________________
 
	% If restriction d=b is imposed, only adjust d.
    if opt.restrictDB
        param = coeffs.db(1);
    else
        param = coeffs.db;
    end

    % Level parameter MuHat.
    if opt.levelParam
        param = [ param reshape(coeffs.muHat,1,p) ];
    end

    % Unrestricted constant.
    if opt.unrConstant
        param = [ param reshape(coeffs.xiHat,1,p)];
    end

	% alphaHat
    if r > 0
        param = [  param reshape( coeffs.alphaHat, 1, p*r ) ];
    end

	% GammaHat
    if k > 0
        param = [  param reshape( coeffs.GammaHat, 1, p*p*k ) ];
    end

end
