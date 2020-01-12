

# Temporary staging area for Matlab code. 

function [ Q, pvQ, LM, pvLM, mvQ, pvMVQ ] <- mv_wntest(x, maxlag, printResults)
# function [ Q, pvQ, LM, pvLM, mvQ, pvMVQ ] <- 
#                       mv_wntest(x, maxlag, printResults)
# Written by Michal Popiel and Morten Nielsen (This version 7.21.2015)
%
# DESCRIPTION: This function performs a multivariate Ljung-Box Q-test for
# 	white noise and univariate Q-tests and LM-tests for white noise on the
# 	columns of x.
# 	The LM test should be consistent for heteroskedastic series, Q-test is not.
%
# Input <- x            (matrix of variables under test, typically model residuals)
#         maxlag       (number of lags for serial correlation tests)
#         printResults (set =1 to print results to screen)
# Output <- Q     (1xp vector of Q statistics for individual series)
%		   pvQ 	 (1xp vector of P-values for Q-test on individual series)
#          LM    (1xp vector of LM statistics for individual series)
%		   pvLM	 (1xp vector of P-values for LM-test on individual series)
#          mvQ   (multivariate Q statistic)
#          pvMVQ (P-value for multivariate Q-statistic using p^2*maxlag df)
%_________________________________________________________________________

p <- size(x,2)

# Create bins for values
pvQ  <- ones(1,p)
pvLM <- ones(1,p)
Q  <- zeros(1,p)
LM <- zeros(1,p)

# Perform univariate Q and LM tests and store the results.
for i <- 1:p
[Q(i),  pvQ(i)] <- Qtest(x(:,i), maxlag)
[LM(i), pvLM(i)] <- LMtest(x(:,i),maxlag)
end

# Perform multivariate Q test.
[mvQ, pvMVQ] <- Qtest(x,maxlag)

# Print output
if printResults
fprintf('\n       White Noise Test Results (lag <- %g)\n', maxlag)
fprintf('---------------------------------------------\n')
fprintf('Variable |       Q  P-val |      LM  P-val  |\n')
fprintf('---------------------------------------------\n')
fprintf('Multivar | %7.3f  %4.3f |     ----  ----  |\n', mvQ, pvMVQ)
for i=1:p
fprintf('Var%g     | %7.3f  %4.3f | %7.3f  %4.3f  |\n',
        i, Q(i), pvQ(i), LM(i), pvLM(i) )
end
fprintf('---------------------------------------------\n')
end
end

function [ LMstat, pv ] <- LMtest(x,q)
# Breusch-Godfrey Lagrange Multiplier test for serial correlation.
T <- size(x,1)
x <- x-mean(x)
y <- x(q+1:end)
z <- x(1:end-q)
for i <- 1:q-1
z <- [x(i+1:end-q+i) z]
end
e <- y
s <- z(:,1:q).*repmat(e,1,q)
sbar <- mean(s)
# The next line bsxfun(@FUNC, A, B) applies the element-by-element binary
# operation FUNC to arrays A and B, with singleton expansion enabled.
s <- bsxfun(@minus, s, sbar)
S <- s'*s/T %'
LMstat <- T*sbar*S^(-1)*sbar' #'
pv <- 1-chi2cdf(LMstat,q)
end


function [ Qstat, pv ] <- Qtest(x, maxlag)
# (Multivariate) Ljung-Box Q-test for serial correlation, see
# 	Luetkepohl (2005, New Introduction to Multiple Time Series Analysis, p. 169).
T <- size(x,1)
p <- size(x,2)

C0=zeros(p)
for t=1:T
C0=C0+x(t,:)'*x(t,:) %'
end
C0=C0./T

C <- zeros(p,p,maxlag) # a3 <- array(seq(12), dim = c(2, 2, 3))
for i=1:maxlag
for t=i+1:T
C(:,:,i)=C(:,:,i)+x(t,:)'*x(t-i,:) %'
end
C(:,:,i)=C(:,:,i)./(T-i) # Note division by (T-i) instead of T.
end

# (Multivariate) Q statistic
Qstat <- 0
for j=1:maxlag
#         Qstat <- Qstat+trace(C(:,:,j)'*inv(C0)*C(:,:,j)*inv(C0)) / (T-j) %'
# The following line is a more efficient calculation than the previous
Qstat <- Qstat+trace( (C(:,:,j)'/C0)*(C(:,:,j)/C0) ) / (T-j) %'
end
Qstat <- Qstat*T*(T+2)
pv <- 1-chi2cdf(Qstat,p*p*maxlag) # P-value is calculated with p^2*maxlag df.
end
