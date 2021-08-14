# Test of knitr

library(knitr)


path <- 'C:/Users/le279259/OneDrive - University of Central Florida/Documents/Research/FCVAR/GitRepo/FCVAR/R_dev/knitr_demo'

setwd(path)


getwd()


# Build the html file with the comments and output.
knitr::spin("FCVAR_demo.R")
