function [ Lbkx ] = Lbk(x, b, k, M)
% 
% function [ Lbkx ] = Lbk(x, b, k, M)
% Lee Morin & Morten Nielsen
% August 21, 2011
%
% Lbk(x, b, k, M) is a lag polynomial in a fractional lag operator.
%
% input = vector or matrix of data x
%       scalar b is the value at which to calculate the fractional difference.
%       scalar k is the number of lags.
%       scalar M is the number of MA lags in the fractional differencing operator.
%           If M is positive, the fractional differencing filter is truncated at lag M.
% 
% output = matrix  [ Lb^1 x, Lb^2 x, ..., Lb^k x]
%       where Lb = 1 - (1-L)^b
% 
% Calls the function xd = FracDiff(x, d, M) 
%______________________________________________________


p = size(x,2);

% Start with Lb^0 x (= x) and overwrite it if k > 0.
Lbkx = x;


% For i = 1, initialize first column of Lbkx = Lb^1 x.
if k > 0
    bx = FracDiff(x, b, M);
    Lbkx = x - bx;
end

if k > 1
    for i = 2:k
    
        Lbkx = [ Lbkx  ( Lbkx(:, p*(i-2)+1 : end ) - FracDiff(Lbkx(:, p*(i-2)+1 : end ), b, M) ) ];

    end
end


% end

