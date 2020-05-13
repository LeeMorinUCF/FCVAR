
##################################################
# 
# Fractionally Cointegrated VAR Model
# Tables for Critical Values and P-values
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
# fracdist_tables.R extracts tables for fracdist 
#   critial values and p-values and organizes them 
#   in a form suitable for use in R. 
# 
# Dependencies:
#   None.
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
# End
##################################################
