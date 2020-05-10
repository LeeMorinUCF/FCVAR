##################################################
#
# Fun with Fortran
#
#
# Calling Fortran subroutines in R.
#
##################################################


##################################################
# Toy example from Wikeversity:
# https://en.wikiversity.org/wiki/R_(programming_language)/...
#   Tutorials/Connecting_Fortran_and_R
##################################################

getwd()

dyn.load("R_dev/Fortran_tests/multiply.dll")


is.loaded("multiply")


# Test an example.
a <- 5
b <- 2
.Fortran("multiply", as.integer(a), as.integer(b), c = integer(1))

# Script to test the link between Fortran and R
star <- function(a,b){
  # res holds the results for the Fortran call
  res <- .Fortran("multiply", as.integer(a), as.integer(b),c = integer(1))
  return(res$c)
}

res_test <- star(a, b)

res_test


##################################################
# Testing fracdiff.f
##################################################


dyn.load("R_dev/Fortran_tests/fracdist_win.dll")

dyn.load("C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR/R_dev/Fortran_tests/fracdist_win.dll")


dyn.load("R_dev/Fortran_tests/fpvalsub.dll")

dyn.load("C:/Users/le279259/Documents/Research/FCVAR/GitRepo/FCVAR/R_dev/Fortran_tests/fpvalsub.dll")



# Test item by item.

# fpval calls olsqc and chicdf

# olsqc calls chol2

# chicdf calls gammp

# gammp calls gser and gcf

# gser calls gammln

# Start from the bottom.


##################################################
# Implementing chol2
##################################################


dyn.load("R_dev/Fortran_tests/chol2.dll")
is.loaded("chol2")

# Arguments: a,m,n,kf2

# A <- matrix(c(2.0, 1.0, 4.0, 1.0), nrow = 2, ncol = 2)

# A_sub <- matrix(c(2.0, 1.0, 2.0, 0.0), nrow = 2, ncol = 2)
A_sub <- matrix(c(2, 1, 2, 0), nrow = 2, ncol = 2)
A_sub %*% t(A_sub)

A <- A_sub %*% t(A_sub)

solve(A)
chol(A)

m <- 2
n <- 2
chol2_out <- .Fortran("chol2", Ainv = as.matrix(A), m = as.integer(m), n = as.integer(n), kf2 = integer(1))

chol2_out$Ainv
# Check.



##################################################
# Implementing olsqc and chol2
##################################################


dyn.load("R_dev/Fortran_tests/olsqc_chol2.dll")
# is.loaded("olsqc_chol2")
is.loaded("olsqc")
is.loaded("chol2")

num_obs <- 100
x_mat <- as.matrix(data.frame(const = rep(1, num_obs),
                              x1 = rnorm(num_obs),
                              x2 = rnorm(num_obs)))
k <- ncol(x_mat)
epsilon <- rnorm(num_obs)
beta <- as.matrix(rep(1, k), nrow = k)
y <- x_mat %*% beta + epsilon

summary(data.frame(cbind(y, x_mat)))

lm_1 <- lm(y ~ x1 + x2, data = data.frame(cbind(y, x_mat)))
summary(lm_1)


# Call the Fortran wrapper.
# Format of subroutine:
# olsqc(nobs,nomax,nvar,nvmax,i1,ssr,beta,
#       & xpy,xpxi,xmat,yvect,resid)

olsqc_out <- .Fortran("olsqc",
                      nobs = as.integer(num_obs),
                      nomax = as.integer(k),
                      nvar = as.integer(k),
                      nvmax = as.integer(k),
                      i1 = as.integer(0),
                      ssr = numeric(1),
                      beta = numeric(k),
                      xpy = matrix(0, nrow = k, ncol = 1),
                      xpxi = matrix(0, nrow = k, ncol = k),
                      xmat = as.matrix(x_mat),
                      yvect = as.matrix(y),
                      resid = matrix(0, nrow = num_obs, ncol = 1))



# Format of subroutine:
# olsqc(nobs,nomax,nvar,nvmax,i1,ssr,beta,
#       & xpy,xpxi,xmat,yvect,resid)


attributes(olsqc_out)

olsqc_out$nobs
olsqc_out$nomax
olsqc_out$nvar
olsqc_out$nvmax
olsqc_out$i1

# Fix these:
olsqc_out$ssr
sqrt(olsqc_out$ssr/(olsqc_out$nobs - k)) # Doesn't match
# Matches the SSresid:
sum(olsqc_out$resid^2)


olsqc_out$beta # Doesn't match
lm_1$coefficients

olsqc_out$xpy
t(x_mat) %*% y # Only the first column matches.
olsqc_out$xpxi
solve(t(x_mat) %*% x_mat)
olsqc_out$xpxi / solve(t(x_mat) %*% x_mat)
solve(t(x_mat) %*% x_mat) / olsqc_out$xpxi
solve(olsqc_out$xpxi)

# The same:
summary(olsqc_out$xmat)
summary(x_mat)
head(olsqc_out$xmat)
head(x_mat)
tail(olsqc_out$xmat)
tail(x_mat)
max(abs(olsqc_out$xmat - x_mat))

# The same:
summary(olsqc_out$yvect)
summary(y)
max(abs(olsqc_out$yvect - y))

# Not the same:
summary(olsqc_out$resid)
summary(lm_1$residuals)

resid_check <- olsqc_out$yvect - olsqc_out$xmat %*% olsqc_out$beta

summary(resid_check)
# Not a match.

# Compare with SSR:
sum(olsqc_out$resid^2)
olsqc_out$ssr
# Match.



# What if its not reaching all observations?
lm_0 <- lm(V1 ~ x1 + x2,
           data = data.frame(cbind(y, x_mat))[seq(k), ])
summary(lm_0)
olsqc_out$beta

##################################################
# End
##################################################


