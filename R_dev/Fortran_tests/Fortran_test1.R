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


# Standard model:
# num_obs <- 100
# x_mat <- as.matrix(data.frame(const = rep(1, num_obs),
#                               x1 = rnorm(num_obs),
#                               x2 = rnorm(num_obs)))
# k <- ncol(x_mat)
# epsilon <- rnorm(num_obs)
# beta <- as.matrix(rep(1, k), nrow = k)
# y <- x_mat %*% beta + epsilon

# Toy example:
num_obs <- 4
x_mat <- as.matrix(data.frame(const = rep(1, num_obs),
                              x1 = sample(num_obs),
                              x2 = sample(num_obs)))
k <- ncol(x_mat)
# Draw the random numbers until round numbers:
epsilon <- sample(num_obs)
# Or just add it to the intercept.

beta <- as.matrix(c(0, rep(1, k-1)), nrow = k)
y <- x_mat %*% beta + epsilon
# Redraw until unique.
y


summary(data.frame(cbind(y, x_mat)))

lm_1 <- lm(y ~ x1 + x2, data = data.frame(cbind(y, x_mat)))
summ_ml_1 <- summary(lm_1)
summ_ml_1

# Call the Fortran wrapper.
# Format of subroutine:
# olsqc(nobs,nomax,nvar,nvmax,i1,ssr,beta,
#       & xpy,xpxi,xmat,yvect,resid)

olsqc_out <- .Fortran("olsqc",
                      nobs = as.integer(num_obs),
                      nomax = as.integer(k),
                      nvar = as.integer(k),
                      nvmax = as.integer(k),
                      # i1 = as.integer(0),
                      i1 = as.integer(1),
                      ssr = numeric(1),
                      beta = numeric(k),
                      xpy = matrix(0, nrow = k, ncol = 1),
                      # xpxi = matrix(0, nrow = k, ncol = k),
                      xpxi = as.matrix(solve(t(x_mat) %*% x_mat)),
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

# Check SSR:
olsqc_out$ssr
# Matches the SSresid:
sum(olsqc_out$resid^2)

# Doesn't match the RSE:
sqrt(olsqc_out$ssr/(olsqc_out$nobs - olsqc_out$nvar))
summ_ml_1$sigma
sqrt(olsqc_out$ssr/(olsqc_out$nobs - olsqc_out$nvar)) /
  summ_ml_1$sigma


olsqc_out$beta # Doesn't match
lm_1$coefficients

# But it does match the calculated xpxi and xpy:
olsqc_out$beta
olsqc_out$xpxi %*% olsqc_out$xpy


# Look at the intermediate calculations:
# xpy
olsqc_out$xpy
t(x_mat) %*% y # Only the first column matches.
t(x_mat) %*% y - olsqc_out$xpy


# xpxi
olsqc_out$xpxi
solve(t(x_mat) %*% x_mat)
olsqc_out$xpxi / solve(t(x_mat) %*% x_mat)
solve(t(x_mat) %*% x_mat) / olsqc_out$xpxi

# Consider the xpx before inversion:
solve(olsqc_out$xpxi)
t(x_mat) %*% x_mat
round(t(x_mat) %*% x_mat - solve(olsqc_out$xpxi), 2)


# The same:
summary(olsqc_out$xmat)
summary(x_mat)
head(olsqc_out$xmat)
head(x_mat)
tail(olsqc_out$xmat)
tail(x_mat)
max(abs(olsqc_out$xmat - x_mat))
summary(abs(olsqc_out$xmat - x_mat))

# The same:
summary(olsqc_out$yvect)
summary(y)
max(abs(olsqc_out$yvect - y))
summary(abs(olsqc_out$yvect - y))

# But these were not changed in the program.


# Not the same:
summary(olsqc_out$resid)
summary(lm_1$residuals)
summary(olsqc_out$resid - lm_1$residuals)


resid_check <- olsqc_out$yvect - olsqc_out$xmat %*% olsqc_out$beta
summary(resid_check)
# Not a match.
# Residuals don't follow from the data and estimates.

# Compare with SSR:
sum(olsqc_out$resid^2)
olsqc_out$ssr
# Match.



# # What if its not reaching all observations?
# lm_0 <- lm(V1 ~ x1 + x2,
#            data = data.frame(cbind(y, x_mat))[seq(k), ])
# summary(lm_0)
# olsqc_out$beta
# Didn't help.

##################################################
# Ridiculously simple example.
##################################################

num_obs <- 4
k <- 3
xpxi <- diag(k)
x_mat <- matrix(0, nrow = num_obs, ncol = k)
x_mat[1, 3] <- 5
y <- matrix(c(7, 0, 0, 0))

t(x_mat) %*% y


olsqc_out <- .Fortran("olsqc",
                      nobs = as.integer(num_obs),
                      nomax = as.integer(k),
                      nvar = as.integer(k),
                      nvmax = as.integer(k),
                      # i1 = as.integer(0),
                      i1 = as.integer(1),
                      ssr = numeric(1),
                      beta = numeric(k),
                      xpy = matrix(0, nrow = k, ncol = 1),
                      # xpxi = matrix(0, nrow = k, ncol = k),
                      xpxi = as.matrix(xpxi),
                      xmat = as.matrix(x_mat),
                      yvect = as.matrix(y),
                      resid = matrix(0, nrow = num_obs, ncol = 1))



# Look at the intermediate calculations:
# Compare xpy
olsqc_out$xpy
t(x_mat) %*% y
# Difference
t(x_mat) %*% y - olsqc_out$xpy



##################################################
# Brute Force on Ridiculously simple example.
##################################################

num_obs <- 4
k <- 3
xpxi <- diag(k)

x_col_num_list <- seq(k)
x_row_num_list <- seq(num_obs)
y_row_num_list <- seq(num_obs)
results <- data.frame(expand.grid(x_col = x_col_num_list,
                                  x_row = x_row_num_list,
                                  y_row = y_row_num_list))
xty_cols <- sprintf('xty_%d', seq(k))
xpy_cols <- sprintf('xpy_%d', seq(k))
diff_cols <- sprintf('diff_%d', seq(k))
results[, xty_cols] <- NA
results[, xpy_cols] <- NA
results[, diff_cols] <- NA


for (combo_num in 1:nrow(results)) {

  x_col <- results[combo_num, 'x_col']
  x_row <- results[combo_num, 'x_row']
  y_row <- results[combo_num, 'y_row']


  # Populate matrices with one element each.
  x_mat <- matrix(0.0, nrow = num_obs, ncol = k)
  x_mat[x_row, x_col] <- combo_num
  # x_mat[x_row, x_col] <- which(matrix(x_mat == combo_num, nrow = num_obs*k))
  x_mat[x_row, x_col] <- which(matrix(x_mat == combo_num,
                                      nrow = num_obs*k, byrow = TRUE))
  y <- matrix(0.0, nrow = num_obs, ncol = 1)
  y[y_row] <- 1

  # Call the Fortran subroutine.
  olsqc_test <- .Fortran("olsqc",
                        nobs = as.integer(num_obs),
                        nomax = as.integer(k),
                        nvar = as.integer(k),
                        nvmax = as.integer(k),
                        # i1 = as.integer(0),
                        i1 = as.integer(1),
                        ssr = numeric(1),
                        beta = numeric(k),
                        xpy = matrix(0, nrow = k, ncol = 1),
                        # xpy = numeric(k),
                        # xpxi = matrix(0, nrow = k, ncol = k),
                        xpxi = as.matrix(xpxi),
                        # xmat = as.matrix(x_mat),
                        # xmat = as.matrix(t(x_mat)),
                        # xmat = as.numeric(x_mat),
                        # xmat = as.numeric(t(x_mat)),
                        # xmat = as.numeric(matrix(x_mat, nrow = num_obs*k, byrow = TRUE)),
                        xmat = as.numeric(matrix(t(x_mat), nrow = num_obs*k, byrow = TRUE)),
                        yvect = as.numeric(y),
                        resid = matrix(0, nrow = num_obs, ncol = 1))



  # Record calculatons and difference.
  results[combo_num, xty_cols] <- t(x_mat) %*% y
  results[combo_num, xpy_cols] <- olsqc_test$xpy
  results[combo_num, diff_cols] <- olsqc_test$xpy - t(x_mat) %*% y


}

# Look for differences.
summary(results)

# Average differences:
colSums(results[, diff_cols])
# On average, they are correct.

# Errors in 16 places.
sum(!(rowSums(results[, diff_cols]) == 0))
sum(!(rowSums(abs(results[, diff_cols])) == 0))

# Any incorrect:
results[!(rowSums(abs(results[, diff_cols])) == 0), ]

# All correct:
# results[(rowSums(abs(results[, diff_cols])) == 0), ]

# Correct a



# Find differences by variable.
var_num <- 1
diff_var <- sprintf('diff_%d', var_num)
results[!(results[, diff_var] == 0), ]
# Column 1 correct.

var_num <- 2
diff_var <- sprintf('diff_%d', var_num)
results[!(results[, diff_var] == 0), ]
# 4 values moved around.

var_num <- 3
diff_var <- sprintf('diff_%d', var_num)
results[!(results[, diff_var] == 0), ]


# Find out where its correct.
var_num <- 1
diff_var <- sprintf('diff_%d', var_num)
results[(results[, diff_var] == 0), ]
sum((results[, diff_var] == 0))
# All entries.

var_num <- 2
diff_var <- sprintf('diff_%d', var_num)
results[(results[, diff_var] == 0), ]
sum(results[, diff_var] == 0)

var_num <- 3
diff_var <- sprintf('diff_%d', var_num)
results[(results[, diff_var] == 0), ]
sum(results[, diff_var] == 0)


##################################################
# End
##################################################


