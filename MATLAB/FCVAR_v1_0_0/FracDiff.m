function [dx] = FracDiff(x, d, M)
%
% function [dx] = FracDiff(x,d,M)
% Lee Morin & Morten Nielsen
% January 31, 2011
%
% FracDiff(x,d,M) is a fractional differencing procedure with optional truncation.
%
% input = vector or matrix of data x
%       scalar d is the value at which to calculate the fractional difference.
%       scalar M is the number of MA lags in the fractional differencing operator.
%           If M is positive, the fractional differencing filter is truncated at lag M.
% 
% output = vector or matrix (1-L)^d x
%
% This program is based on code borrowed from Katsumi Shimotsu, September 2003.
%______________________________________________________


T = size(x, 1);

if M > 0
    k = (1 : M-1)';
else
    k = (1 : T-1)';
end

b = (k-d-1)./k;
b = [1; cumprod(b)];
dx = filter(b,1,x);



% end