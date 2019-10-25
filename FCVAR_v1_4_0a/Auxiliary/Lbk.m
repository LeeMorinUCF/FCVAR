
function [ Lbkx ] = Lbk(x, b, k)
% function [ Lbkx ] = Lbk(x, b, k)
% Written by Michal Popiel and Morten Nielsen (This version 10.22.2014)
% Based on Lee Morin & Morten Nielsen (May 24, 2013)
%
% DESCRIPTION: Lbk(x, b, k) is a lag polynomial in the fractional lag operator.
%
% Input = x (vector or matrix of data)
%         b (scalar value at which to calculate the fractional lag)
%         k (number of lags)
% 
% Output = matrix  [ Lb^1 x, Lb^2 x, ..., Lb^k x] where Lb = 1 - (1-L)^b. The output 
%		matrix has the same number of rows as x but k times as many columns.
% 
% Calls the function FracDiff(x, d) 
%______________________________________________________


    p = size(x,2);
    
    % Initialize output matrix.
    Lbkx = [];

    % For i = 1, set first column of Lbkx = Lb^1 x.
    if k > 0
        bx = FracDiff(x, b);
        Lbkx = x - bx;
    end

    if k > 1
        for i = 2:k   
            Lbkx = [ Lbkx  ( Lbkx(:, p*(i-2)+1 : end ) - FracDiff(Lbkx(:, p*(i-2)+1 : end ), b) ) ];

        end
    end
end