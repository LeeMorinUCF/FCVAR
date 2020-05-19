
##################################################
# 
# Fractionally Cointegrated VAR Model
# Test Cases for Critical Values and P-values
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
# fracdist_tables.R creates tables for fracdist 
#   critial values and p-values and organizes them 
#   in a form suitable for evaluation in Fortran. 
# 
# Dependencies:
#   None.
# 
# TODO: Revise output file headings for csv file format. 
# 
##################################################


##################################################
# Preparing the Workspace
##################################################

# Clear workspace.
rm(list=ls(all=TRUE))

# Set working directory.
wd_path <- '~/Research/FCVAR'

# Set data directory.
data_dir <- '~/Research/FCVAR/fracdist/mn-files'

setwd(wd_path)


##################################################
# Load Packages
##################################################

# Inspired by a subroutine in the Fortran fracdist code 
# subroutine readdata(iq,iscon,probs,bbb,xndf)
# c This routine reads data from whichever input file is appropriate
# c for the no-constant case.
# c If iscon=0, files are frmapp01.txt through frmapp12.txt
# c If iscon.ne.0, files are frcapp01.txt through frcapp12.txt




# This function assembles a table of probabilities 
# and quantiles for all values of the fractional 
# integration parameter for a particular rank and 
# soecification of constant. 
get_fracdist_tab <- function(iq, iscon, dir_name) {
  
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
  
  
  return(frtab)
}



frtab <- get_fracdist_tab(iq = 1, iscon = 0, dir_name = data_dir)

summary(frtab)
head(frtab)
frtab[217:223, ]
tail(frtab)


##################################################
# Loop over tables and save them in a compressed
# format suitable for R packages. 
##################################################


# Set directory for output.
out_dir <- '~/Research/FCVAR/fracdist/R-mn-files'

for (iq in seq(12)) {
  
  for (iscon in c(0, 1)) {
    
    # Read text file version of table from Fortran package.
    frtab <- get_fracdist_tab(iq = 1, iscon = 0, dir_name = data_dir)
    
    
    # Determine corresponding file name for output. 
    
    # Depends on whether or not there is a constant term.
    if (iscon == 0) {
      dfirst = 'frmapp'
    } else {
      dfirst = 'frcapp'
    }
    
    # Depends on the (difference in) cointerating ranks. 
    dq <- sprintf('00%d', iq)
    dq <- substr(dq, nchar(dq) - 1, nchar(dq))
    
    
    out_file_name <- sprintf('%s/%s%s.RData', out_dir, dfirst, dq)
    
    
    # Save the file in a compressed format suitable for R. 
    save(frtab, file = out_file_name)
    
  }
  
}


# After initial save, test the optimal compression format. 
tools::checkRdaFiles(out_dir)
# Looks good. 54K each file, 1.24M total. 

# Inspect one file.
check_frtab <- get(load(file = out_file_name))

class(check_frtab)
head(check_frtab)
tail(check_frtab)
summary(check_frtab)

# See revised version of get_fracdist_tab function:
# get_fracdist_tab(iq, iscon, dir_name, file_ext = 'txt')



##################################################
# Generate Test Cases
##################################################

# Set lists of input variables.
iscon_list <- c(0, 1)
iq_list <- seq(12)
clevel_list <- c(0.10, 0.05, 0.01)
num_bb <- 10
num_stat <- 10
set.seed(42)
bb_list <- runif(num_bb, min = 0.50, max = 2.0)
stat_inv_list <- runif(num_stat, min = 0.80, max = 1.0)


# Create two tables, one for each function.


# Test p-values for a variety of values of the statistic.
test_fpval <- expand.grid(bb = bb_list, 
                          iscon = iscon_list, 
                          stat = rep(NA, num_stat), 
                          iq = iq_list)
# Draw test statistic values from the corresponding chi-squared distribution. 
iq_length <- 2*num_bb*num_stat
for (iq_num in 1:length(iq_list)) {
  
  iq <- iq_list[iq_num]
  row_sel <- seq((iq_num - 1)*iq_length + 1, iq_num*iq_length)
  test_fpval[row_sel, 'stat'] <- qchisq(p = rep(stat_inv_list, 2*num_bb), 
                                        df = iq^2)
  
}
# Reorder the columns. 
test_fpval <- test_fpval[, c('iscon', 'iq', 'bb', 'stat')]
summary(test_fpval)
head(test_fpval)
tail(test_fpval)


# Test critical values for conventional significance levels. 
test_fcval <- expand.grid(bb = bb_list, 
                          iscon = iscon_list, 
                          clevel = clevel_list, 
                          iq = iq_list)
test_fcval <- test_fcval[, c('iscon', 'iq', 'bb', 'clevel')]
# test_fcval <- test_fcval[order(test_fcval$), ]
summary(test_fcval)
head(test_fcval)
tail(test_fcval)


# Save the files in fixed-width format. 
# Yes, I know this is slow but I want to control the formatting
# the way I would read it in Fortran. 
out_file_name <- 'test_fpval.txt'
out_file_name <- sprintf('%s/%s', out_dir, out_file_name)
cat(sprintf('%s\n', paste(colnames(test_fpval), collapse = ' ')), 
    file = out_file_name)
for (line in 1:nrow(test_fpval)) {
  
  cat(sprintf('%d ', test_fpval[line, 'iscon']), 
      file = out_file_name, append = TRUE)
  iq <- test_fpval[line, 'iq']
  if (iq >= 10) {
    cat(sprintf('%d ', iq), 
        file = out_file_name, append = TRUE)
  } else {
    cat(sprintf(' %d ', iq), 
        file = out_file_name, append = TRUE)
  }
  cat(sprintf('%5.3f ', test_fpval[line, 'bb']), 
      file = out_file_name, append = TRUE)
  cat(sprintf('%8.4f\n', test_fpval[line, 'stat']), 
      file = out_file_name, append = TRUE)
  
}


# Save the files in fixed-width format. 
out_file_name <- 'test_fcval.txt'
out_file_name <- sprintf('%s/%s', out_dir, out_file_name)
cat(sprintf('%s\n', paste(colnames(test_fcval), collapse = ' ')), 
    file = out_file_name)
for (line in 1:nrow(test_fcval)) {
  
  cat(sprintf('%d ', test_fcval[line, 'iscon']), 
      file = out_file_name, append = TRUE)
  iq <- test_fcval[line, 'iq']
  if (iq >= 10) {
    cat(sprintf('%d ', iq), 
        file = out_file_name, append = TRUE)
  } else {
    cat(sprintf(' %d ', iq), 
        file = out_file_name, append = TRUE)
  }
  cat(sprintf('%5.3f ', test_fcval[line, 'bb']), 
      file = out_file_name, append = TRUE)
  cat(sprintf('%4.2f\n', test_fcval[line, 'clevel']), 
      file = out_file_name, append = TRUE)
  
}


##################################################
# End
##################################################
