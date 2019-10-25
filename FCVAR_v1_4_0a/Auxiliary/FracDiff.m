function [dx] = FracDiff(x, d)
% function [dx] = FracDiff(x,d)
% Andreas Noack Jensen & Morten Nielsen
% May 24, 2013
%
% FracDiff(x,d) is a fractional differencing procedure based on the
% 	fast fractional difference algorithm of Jensen & Nielsen (2014, JTSA).
%
% input = x (vector or matrix of data)
%         d (scalar value at which to calculate the fractional difference)
% 
% output = vector or matrix (1-L)^d x of same dimensions as x.
%______________________________________________________


    [T, p] = size(x);
    k = (1 : T-1)';

    % NEXTPOW2(N) returns the first P such that 2.^P >= abs(N).  It is
    %     often useful for finding the nearest power of two sequence
    %     length for FFT operations.
    NFFT = 2^nextpow2(2*T-1);

    % Array operation of the index of the series without the last element minus
    % the order of integration+1, divided by that same series
    b = (k-d-1)./k;

    % cumulative product of that series modified in previous line
    b = [1; cumprod(b)];

    % IFFT(X) is the inverse discrete Fourier transform of X.
    %     IFFT(..., 'symmetric') causes IFFT to treat X as conjugate symmetric
    %     along the active dimension.  This option is useful when X is not exactly
    %     conjugate symmetric merely because of round-off error. 
    % REPMAT Replicate and tile an array.
    %     B = repmat(A,M,N) creates a large matrix B consisting of an M-by-N
    %     tiling of copies of A. The size of B is [size(A,1)*M, size(A,2)*N].
    %     The statement repmat(A,N) creates an N-by-N tiling.
    dx = ifft(repmat(fft(b, NFFT), 1, p).*fft(x, NFFT), 'symmetric');

    dx = dx(1:T,:);
end