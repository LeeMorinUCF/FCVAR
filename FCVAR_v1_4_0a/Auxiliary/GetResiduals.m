function [ epsilon ] = GetResiduals(x, k, r, coeffs, opt)
% function [ epsilon ] = GetResiduals(x, k, r, coeffs, opt)
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
% Based on Lee Morin & Morten Nielsen (August 22, 2011)
% 
% DESCRIPTION: This function calculates the model residuals.
%
% Input = x (matrix of variables to be included in the system)
%         k (number of lags)
%         r (number of cointegrating vectors)
%         coeffs (Matlab structure of coefficients)
%         opt (object containing the estimation options)
% Output = epsilon (matrix of residuals from model estimation evaluated at
%                     the parameter estimates specified in coeffs) 
%_________________________________________________________________________


    % If level parameter is included, the data must be shifted before
    %   calculating the residuals:
    if opt.levelParam
        T = size(x,1);
        y = x - ones(T,1)*coeffs.muHat;
    else
        y = x;
    end

    % --- Transform data --- %
    [ Z0, Z1, Z2, Z3 ] = TransformData(y, k, coeffs.db, opt);

    % --- Calculate residuals --- %
    epsilon = Z0;

    if r > 0
        epsilon = epsilon - Z1*[ coeffs.betaHat; coeffs.rhoHat ]*coeffs.alphaHat';
    end

    if k > 0
        epsilon = epsilon - Z2*coeffs.GammaHat';
    end
    
    if opt.unrConstant
        epsilon = epsilon - Z3*coeffs.xiHat';
    end
        
end

