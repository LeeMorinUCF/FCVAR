

# Temporary staging area for Matlab code. 


function [ like ] <- FCVARlike(x, params, k, r, opt)
# function [ like ] <- FCVARlike(x, params, k, r, opt)
# Written by Michal Popiel and Morten Nielsen (This version 11.10.2014)
# 
# DESCRIPTION: This function adjusts the variables with the level parameter,
# 	if present, and returns the log-likelihood given d,b.
%
# Input <- x   (matrix of variables to be included in the system)
#         params (a vector of parameters d,b, and mu (if option selected))
#         k   (number of lags)
#         r   (number of cointegrating vectors)
#         opt (object containing the estimation options)
# Output <- like (concentrated log-likelihood evaluated at given parameters)
%_________________________________________________________________________
global estimatesTEMP

T <- size(x,1)
p <- size(x,2)


# If there is one linear restriction imposed, optimization is over phi,
#  so translate from phi to (d,b), otherwise, take parameters as
#  given.
if(size(opt$R_psi,1)==1)
  H <- null(opt$R_psi)
h <- opt$R_psi'*inv(opt$R_psi*opt$R_psi')*opt$r_psi       
db <- H*params(1)+h
d <- db(1)
b <- db(2)
if (opt$levelParam)
  mu <- params(2:end)
end
else
  d <- params(1)
b <- params(2)    
if (opt$levelParam)
  mu <- params(3:end)
end
end

# Modify the data by a mu level.
if (opt$levelParam)
  y <- x - ones(T,p)*diag(mu)
else
  y <- x
end

# If d=b is not imposed via opt$restrictDB, but d>=b is required in
%	opt$constrained, then impose the inequality here.
if(opt$constrained)
  b <- min(b,d)
end

# Obtain concentrated parameter estimates$
[ estimates ] <- GetParams(y, k, r, [d b], opt)

# Storing the estimates in a global structure here allows us to skip a
#   call to GetParams after optimization to recover the coefficient
#   estimates
estimatesTEMP <- estimates
# If level parameter is present, they are the last p parameters in the
#   params vector
if (opt$levelParam)
  estimatesTEMP$muHat <- mu
else
  estimatesTEMP$muHat <- NULL
end

# Calculate value of likelihood function.
T <- size(y,1) - opt$N
p <- size(y,2)
like <- - T*p/2*( log(2*pi) + 1)  - T/2*log(det(estimates$OmegaHat))

end
