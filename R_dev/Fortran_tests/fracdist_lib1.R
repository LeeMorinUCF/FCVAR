
##################################################
# 
# Fractionally Cointegrated VAR Model
# Functions for Critical Values and P-values
# 
# Lealand Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# May 12, 2020
# 
##################################################
# 
# fracdist_lib.R calculates critial values and p-values 
#   for rank tests in the FCVAR model. 
# 
# Dependencies:
#   None.
# 
##################################################


##################################################
# Preparing the Workspace
##################################################

# # Clear workspace.
# rm(list=ls(all=TRUE))
# 
# # Set working directory.
# wd_path <- '~/Research/FCVAR'
# 
# # Set data directory.
# data_dir <- '~/Research/FCVAR/fracdist/mn-files'
# 
# setwd(wd_path)


##################################################
# Load Packages
##################################################

# # Generate filename.
# iq <- 1
# iscon <- 0
# bval <- 0.51
# 
# # subroutine readdata(iq,iscon,probs,bbb,xndf)
# # c This routine reads data from whichever input file is appropriate
# # c for the no-constant case.
# # c If iscon=0, files are frmapp01.txt through frmapp12.txt
# # c If iscon.ne.0, files are frcapp01.txt through frcapp12.txt
# 
# if (iscon == 0) {
#   dfirst = 'frmapp'
# } else {
#   dfirst = 'frcapp'
# }
# 
# dq <- sprintf('00%d', iq)
# dq <- substr(dq, nchar(dq) - 1, nchar(dq))
# 
# 
# in_file_name <- sprintf('%s/%s%s.txt', data_dir, dfirst, dq)
# 
# 
# # frtab <- read.table(in_file_name)
# # read.fortran()
# 
# 
# # Let's go old school
# # do ib=1,31
# # read(2,200) bbb(ib)
# # 200    format(4x,f5.2)
# # do j=1,221
# # read(2,201) probs(j), xndf(j,ib)
# # 201      format(f6.4,3x,f16.12)
# # end do
# # end do
# 
# frtab <- data.frame(bbb = numeric(221*31), 
#                     probs = numeric(221*31), 
#                     xndf = numeric(221*31))
# 
# frtab_lines <- readLines(con = in_file_name)
# lines_read <- 0
# 
# # ib <- 1
# # ib <- 2
# for (ib in 1:31) {
#   
#   print(sprintf('Reading for ib = %d.', ib))
#   
#   # Read the value of b.
#   # bbb_str <- readLines(con = in_file_name, n = 1)
#   bbb_str <- frtab_lines[lines_read + 1]
#   lines_read <- lines_read + 1
#   bbb <- as.numeric(substr(bbb_str, 6, 9))
#   
#   # Read the probabilities and quantiles.
#   frtab_sub <- data.frame(bbb = rep(bbb, 221), 
#                           probs = numeric(221), 
#                           xndf = numeric(221))
#   
#   frtab_lines_sub <- frtab_lines[seq(lines_read + 1, lines_read + 221)]
#   lines_read <- lines_read + 221
#   frtab_sub[, 'probs'] <- as.numeric(substr(frtab_lines_sub, 1, 6))
#   frtab_sub[, 'xndf'] <- as.numeric(substr(frtab_lines_sub, 9, 25))
#   
#   frtab[seq((ib - 1)*221 + 1, ib*221), ] <- frtab_sub
#   
#   
#   print(sprintf('lines_read = %d.', lines_read))
#   
# }
# 
# 
# 
# 
# summary(frtab)
# head(frtab)
# frtab[217:223, ]
# tail(frtab)


# Create a function to assemble a table of probabilities 
# and quantiles for all values of the fractional 
# integration parameter for a particular rank and 
# soecification of constant. 
get_fracdist_tab <- function(iq, iscon, dir_name, file_ext = 'txt') {
  
  # Determine required file name. 
  
  # Depends on whether or not there is a constant term.
  if (iscon == 0) {
    dfirst = 'frmapp'
  } else {
    dfirst = 'frcapp'
  }
  
  # Depends on the (difference in) cointerating ranks. 
  dq <- sprintf('00%d', iq)
  dq <- substr(dq, nchar(dq) - 1, nchar(dq))
  
  
  if (file_ext == 'txt') {
    
    # Load original tables from Fortran version.
    
    in_file_name <- sprintf('%s/%s%s.txt', dir_name, dfirst, dq)
    
    # Initialize matrix and open connection.
    frtab <- data.frame(bbb = numeric(221*31), 
                        probs = numeric(221*31), 
                        xndf = numeric(221*31))
    
    frtab_lines <- readLines(con = in_file_name)
    lines_read <- 0
    
    # Loop on fractional integration parameter
    # to collect table from separate pieces.
    for (ib in 1:31) {
      
      # Read the value of b.
      bbb_str <- frtab_lines[lines_read + 1]
      lines_read <- lines_read + 1
      bbb <- as.numeric(substr(bbb_str, 6, 9))
      
      # Read the probabilities and quantiles.
      frtab_sub <- data.frame(bbb = rep(bbb, 221), 
                              probs = numeric(221), 
                              xndf = numeric(221))
      
      frtab_lines_sub <- frtab_lines[seq(lines_read + 1, lines_read + 221)]
      lines_read <- lines_read + 221
      frtab_sub[, 'probs'] <- as.numeric(substr(frtab_lines_sub, 1, 6))
      frtab_sub[, 'xndf'] <- as.numeric(substr(frtab_lines_sub, 9, 25))
      
      # Append to the full table. 
      frtab[seq((ib - 1)*221 + 1, ib*221), ] <- frtab_sub
      
    }
    
  } else if (file_ext == 'RData') {
    
    # Load compressed files for R.
    in_file_name <- sprintf('%s/%s%s.RData', dir_name, dfirst, dq)
    check_frtab <- get(load(file = in_file_name))
    
  } else {
    stop('File extension not supported.')
  }
  
  
  return(frtab)
}



# frtab <- get_fracdist_tab(iq = 1, iscon = 0, dir_name = data_dir)
# 
# summary(frtab)
# head(frtab)
# frtab[217:223, ]
# tail(frtab)

##################################################
# Interpolate critical values local to chosen b.
##################################################

# blocal <- function(nb, bb, estcrit, bval)

blocal <- function(nb, bb, estcrit, bval) {
  
  if (bb < 0.51 | bb > 2.0) {
    stop(sprintf('Specified value of b = %f is out of range 0.51 to 2.0.', bb))
  }
  
  # Weights on neighboring quantiles are based on trapezoidal kernel. 
  weight <- 1.0 - 5.0*abs(bval - bb)
  weight[weight < 0] <- 0
  # Includes at most 9 observations. 
  
  # Determine the quantiles to be used for interpolation. 
  # jbot will be index of lowest b that gets positive weight. 
  jbot <- which(weight > 0)[1]
  nobs <- sum(weight > 0)
  weight_sel <- seq(jbot, jbot + nobs - 1)
  
  # Create a dataset for interpolation by regression.
  yx_mat <- data.frame(y = numeric(9), 
                       x1 = numeric(9), 
                       x2 = numeric(9), 
                       x3 = numeric(9))
  
  yx_mat[1:nobs, 'y'] <- weight[weight_sel]*estcrit[weight_sel]
  yx_mat[1:nobs, 'x1'] <- weight[weight_sel]
  yx_mat[1:nobs, 'x2'] <- weight[weight_sel]*bval[weight_sel]
  yx_mat[1:nobs, 'x3'] <- weight[weight_sel]*bval[weight_sel]^2
  
  lm_olsqc <- lm(formula = y ~ x1 + x2 + x3 - 1, data = yx_mat)
  
  bfit <- sum(lm_olsqc$coefficients*bb^seq(0,2))
  
  return(bfit)
}


##################################################
# Intermediate calculations for p-values
##################################################


# real*8 xndf(221,31), bval(31), probs(221), ginv(221)
# real*8 bedf(221), estcrit(31)
# do i=1,np
#   do j=1,nb
#     estcrit(j) = xndf(i,j)
#   end do
#   call blocal(nb,bb,estcrit,bval,bedf(i))
# end do

# nb <- 31
# np <- 221
# 
# bval <- unique(frtab[, 'bbb'])
# bval <- bval[order(bval)]
# 
# probs <- unique(frtab[, 'probs'])
# probs <- probs[order(probs)]
# 
# bb <- 0.73
# 
# bedf <- rep(NA, np)
# 
# # ib <- 1
# for (i in 1:np) {
#   
#   # b_i <- bval[ib]
#   # estcrit <- frtab[frtab[, 'bbb'] == b_i, 'xndf']
#   
#   prob_i <- probs[i]
#   estcrit <- frtab[frtab[, 'probs'] == prob_i, 'xndf']
#   
#   bedf[i] <- blocal(nb, bb, estcrit, bval)
#   
# }
# 
# # Obtain inverse CDF of the Chi-squared distribution.
# # gcinv(iq,np,probs,ginv)
# ginv <- qchisq(probs, df = iq^2) 


# Next, call the fuction to calculate p-values:
# fpval(np = 9, iq, stat, probs, bedf, ginv, pval)

##################################################
# Calculate p-values
##################################################

# stat <- 3.84

# c This routine calculates P values.
# c
# c stat is test statistic.
# c npts is number of points for local approximation (probably 9).
# c bedf contains quantiles of numerical distribution for specified.
# c value of b or values of b and d.
# c ginv contains quantiles of approximating chi-squared distribution.
# c
fpval <- function(npts = 9, iq, stat, probs, bedf, ginv) {
  
  # nomax <- 25
  # nvmax <- 3
  ndf <- iq**2
  
  # Deal with extreme cases.
  btiny <- 0.5*bedf[1]
  bbig <- 2.0*bedf[221]
  if (stat < btiny) {
    pval <- 1.0
    return(pval)
  }
  if (stat > bbig) {
    pval <- 0.0
    return(pval)
  }
  
  
  # Find critical value closest to test statistic.
  diff <- abs(stat - bedf)
  imin <- which.min(diff)[1]
  diffm <- diff[imin]
  
  # nph <- npts/2 # Division by integer does floor. 
  nph <- floor(npts/2)
  nptop <- 221 - nph
  
  
  # Create a dataset for interpolation by regression.
  yx_mat <- data.frame(y = numeric(npts), 
                       x1 = numeric(npts), 
                       x2 = numeric(npts), 
                       x3 = numeric(npts))
  # The form depends on the proximity to the endpoints. 
  
  if (imin > nph & imin < nptop) {
    # imin is not too close to the end. Use np points around stat.
    
    
    # Populate the dataset for interpolation by regression.
    ic <- imin - nph - 1 + seq(1, npts)
    
    
    # print('imin = ')
    # print(imin)
    # print('nph = ')
    # print(nph)
    # print('npts = ')
    # print(npts)
    # 
    # np1 <- length(ic)
    # 
    # print('np1 = ')
    # print(np1)
    # print('ic = ')
    # print(ic)
    # print('ginv[ic] = ')
    # print(ginv[ic])
    
    yx_mat[, 'y'] <- ginv[ic]
    yx_mat[, 'x1'] <- 1.0
    yx_mat[, 'x2'] <- bedf[ic]
    yx_mat[, 'x3'] <- bedf[ic]^2
    
  } else {
    # imin is close to one of the ends. Use points from imin +/- nph to end.
    
    if (imin < nph) {
      
      np1 <- imin + nph
      np1 <- max(np1, 5)
      
      # Populate the dataset for interpolation by regression.
      ic <- seq(1, np1)
      yx_mat[1:np1, 'y'] <- ginv[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- bedf[ic]
      yx_mat[1:np1, 'x3'] <- bedf[ic]^2
      
    } else { 
      
      np1 <- 222 - imin + nph
      np1 <- max(np1, 5)
      
      # Populate the dataset for interpolation by regression.
      ic <- 222 - seq(1, np1)
      yx_mat[1:np1, 'y'] <- ginv[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- bedf[ic]
      yx_mat[1:np1, 'x3'] <- bedf[ic]^2
      
    }
    
  }
  
  # Run regression and obtain p-value. 
  lm_olsqc <- lm(formula = y ~ x1 + x2 + x3 - 1, data = yx_mat)
  
  crfit <- sum(lm_olsqc$coefficients*stat^seq(0,2))
  crfit <- max(crfit, 10^(-6))
  
  pval <- pchisq(crfit, df = ndf) 
  pval <- 1.0 - pval
  
  return(pval)
}

# pval <- fpval(npts = 9, iq, stat, probs, bedf, ginv)



##################################################
# Main function for calculating p-values
# Includes preliminary calculations. 
##################################################

fracdist_pvalues <- function(iq, iscon, dir_name, bb, stat) {
  
  # Obtain relevant table of statistics.
  frtab <- get_fracdist_tab(iq, iscon, dir_name)
  bval <- unique(frtab[, 'bbb'])
  bval <- bval[order(bval)]
  probs <- unique(frtab[, 'probs'])
  probs <- probs[order(probs)]
  
  # Calculate the approximate CDF for this particular value of b.
  nb <- 31
  np <- 221
  bedf <- rep(NA, np)
  for (i in 1:np) {
    
    prob_i <- probs[i]
    estcrit <- frtab[frtab[, 'probs'] == prob_i, 'xndf']
    
    bedf[i] <- blocal(nb, bb, estcrit, bval)
    
  }
  
  # Obtain inverse CDF of the Chi-squared distribution.
  ginv <- qchisq(probs, df = iq^2) 
  # This is the dependent variable in the interpolation regressions. 
  
  # Calculate p-values. 
  pval <- fpval(npts = 9, iq, stat, probs, bedf, ginv)
  
  return(pval)
}


# pval <- fracdist_pvalues(iq = 1, iscon = 0, dir_name = data_dir, 
#                          bb = 0.73, stat = 3.84)
# print(pval)


##################################################
# Might as well finish the job
# Critical Values, in case someone needs them
##################################################

# c This routine calculates critical values.
# c
# c clevel is level for test.
# c npts is number of points for local approximation (probably 9).
# c bedf contains quantiles of numerical distribution for specified
# c value of b or values of b and d.
# c ginv contains quantiles of approximating chi-squared distribution.
# c
fpcrit <- function(npts = 9, iq, clevel, probs, bedf, ginv) {
  
  # nomax <- 25
  # nvmax <- 3
  ndf <- iq**2
  
  # Handle extreme cases.
  ptiny <- 0.0001
  pbig <- 0.9999
  if (clevel < ptiny) {
    ccrit <- bedf[221]
    return(ccrit)
  }
  if (clevel > pbig) {
    ccrit <- bedf[221]
    return(ccrit)
  }
  
  
  # Find probability closest to test level. 
  cquant <- 1.0 - clevel
  diff <- abs(cquant - probs)
  imin <- which.min(diff)[1]
  diffm <- diff[imin]
  
  gcq <- qchisq(cquant, df = ndf) 
  
  nph <- npts/2
  nptop <- 221 - nph
  
  
  # Create a dataset for interpolation by regression.
  yx_mat <- data.frame(y = numeric(npts), 
                       x1 = numeric(npts), 
                       x2 = numeric(npts), 
                       x3 = numeric(npts))
  # The form depends on the proximity to the endpoints. 
  
  if (imin > nph & imin < nptop) {
    # imin is not too close to the end. Use np points around stat.
    
    # Populate the dataset for interpolation by regression.
    ic <- imin - nph - 1 + seq(1, npts)
    yx_mat[, 'y'] <- bedf[ic]
    yx_mat[, 'x1'] <- 1.0
    yx_mat[, 'x2'] <- ginv[ic]
    yx_mat[, 'x3'] <- ginv[ic]^2
    
  } else {
    # imin is close to one of the ends. Use points from imin +/- nph to end.
    
    if (imin < nph) {
      
      np1 <- imin + nph
      np1 <- max(np1, 5)
      
      # Populate the dataset for interpolation by regression.
      ic <- seq(1, np1)
      yx_mat[1:np1, 'y'] <- bedf[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- ginv[ic]
      yx_mat[1:np1, 'x3'] <- ginv[ic]^2
      
    } else { 
      
      np1 <- 222 - imin + nph
      np1 <- max(np1, 5)
      
      # Populate the dataset for interpolation by regression.
      ic <- 222 - seq(1, np1)
      yx_mat[1:np1, 'y'] <- bedf[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- ginv[ic]
      yx_mat[1:np1, 'x3'] <- ginv[ic]^2
      
    }
    
  }
  
  # Run regression and obtain critical value. 
  lm_olsqc <- lm(formula = y ~ x1 + x2 + x3 - 1, data = yx_mat)
  
  ccrit <- sum(lm_olsqc$coefficients*gcq^seq(0,2))
  
  return(ccrit)
}



# ccrit <- fpcrit(npts = 9, iq, clevel = 0.05, probs, bedf, ginv)
# print(ccrit)


# for (clevel_i in c(0.01, 0.05, 0.10)) {
#   ccrit <- fpcrit(npts = 9, iq, clevel = clevel_i, probs, bedf, ginv)
#   print(ccrit)
#   
# }



##################################################
# Main function for calculating either critical 
# values or p-values
# Includes preliminary calculations. 
##################################################

fracdist_values <- function(iq, iscon, dir_name, bb, stat, 
                            ipc = TRUE, clevel = c(0.01, 0.05, 0.10)) {
  
  # Obtain relevant table of statistics.
  frtab <- get_fracdist_tab(iq, iscon, dir_name)
  bval <- unique(frtab[, 'bbb'])
  bval <- bval[order(bval)]
  probs <- unique(frtab[, 'probs'])
  probs <- probs[order(probs)]
  
  # Calculate the approximate CDF for this particular value of b.
  nb <- 31
  np <- 221
  bedf <- rep(NA, np)
  for (i in 1:np) {
    
    prob_i <- probs[i]
    estcrit <- frtab[frtab[, 'probs'] == prob_i, 'xndf']
    
    bedf[i] <- blocal(nb, bb, estcrit, bval)
    
  }
  
  # Obtain inverse CDF of the Chi-squared distribution.
  ginv <- qchisq(probs, df = iq^2) 
  # This is the dependent variable in the interpolation regressions. 
  
  # Perform required calculation.
  if (ipc == TRUE) {
    
    # Calculate p-values. 
    outval <- fpval(npts = 9, iq, stat, probs, bedf, ginv)
    
  } else {
    
    # Calculate critical values. 
    
    # Might specify a list of critical values. 
    outval <- rep(NA, length(clevel))
    for (i in 1:length(clevel)) {
      clevel_i <- clevel[i]
      outval[i] <- fpcrit(npts = 9, iq, clevel = clevel_i, probs, bedf, ginv)
      
    }
    
  }
  
  return(outval)
}


##################################################
# Testing
##################################################



# Testing critical values.


# With precalculated values from given bb and iq. 
# ccrit <- fpcrit(npts = 9, iq, clevel = 0.05, probs, bedf, ginv)
# print(ccrit)
# 
# for (clevel_i in c(0.01, 0.05, 0.10)) {
#   ccrit <- fpcrit(npts = 9, iq, clevel = clevel_i, probs, bedf, ginv)
#   print(ccrit)
# }

# fracdist_out <- fracdist_values(iq, iscon, dir_name, bb, stat, 
#                                 ipc = FALSE, clevel = 0.05)
# print(fracdist_out)
# 
# fracdist_out <- fracdist_values(iq, iscon, dir_name, bb, stat, 
#                                 ipc = FALSE)
# print(fracdist_out)


# Testing p-values.

# # Using function specifically for p-values:
# pval <- fracdist_pvalues(iq = 1, iscon = 0, dir_name = data_dir, 
#                          bb = 0.73, stat = 3.84)
# print(pval)
# 
# pval <- fracdist_pvalues(iq = 12, iscon = 0, dir_name = data_dir, 
#                          bb = 0.73, stat = 84)
# print(pval)



# fracdist_out <- fracdist_values(iq = 1, iscon = 0, dir_name = data_dir, 
#                                 bb = 0.73, stat = 3.84)
# print(fracdist_out)
# 
# fracdist_out <- fracdist_values(iq = 12, iscon = 0, dir_name = data_dir, 
#                                 bb = 0.73, stat = 3.84)
# print(fracdist_out)
# 
# fracdist_out <- fracdist_values(iq = 12, iscon = 0, dir_name = data_dir, 
#                                 bb = 0.73, stat = 84)
# print(fracdist_out)
# # Fail: Gives p-value for 3.84. Should be 1.00. 

# fracdist_out <- fracdist_values(iq = 12, iscon = 0, dir_name = data_dir, 
#                                 bb = 0.73, stat = 184)
# print(fracdist_out)





##################################################
# End
##################################################
