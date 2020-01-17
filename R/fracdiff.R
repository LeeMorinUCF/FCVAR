# fracdiff <- function(x, d)
# Andreas Noack Jensen & Morten Ørregaard Nielsen
# May 10, 2013
#
# fracdiff(x, d) is a fractional differencing procedure based on the fast fractional
#	differencing algorithm in
#
# Jensen and Nielsen (2013). A fast fractional differencing algorithm.
#	QED working paper 1307, Queen's University.
#
# input = vector of data x
#       scalar d is the value at which to calculate the fractional difference.
# 
# output = vector (1-L)^d x of same dimension as x.

fracdiff_JN <- function(x, d){
  
  ##################################################
  # In the fracdiff package they demean x first.
  # x <- x - mean(x)
  # With the above line added, the results are the same.
  ##################################################
  
  iT <- length(x)
  np2 <- nextn(2*iT - 1, 2)
  k <- 1:(iT-1)
  b <- c(1, cumprod((k - d - 1)/k))
  dx <- fft(fft(c(b, rep(0, np2 - iT))) * fft(c(x, rep(0, np2 - iT))), inverse = T) / np2;
  return(Re(dx[1:iT]))
}