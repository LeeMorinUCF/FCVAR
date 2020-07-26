

# Update R version.
install.packages('installr')
library(installr)

updateR()


# Testing functions to locate local optima.

install.packages('ggpmisc')
library(ggpmisc)

test_1 <- c(1,1,2,3,2,1,1,5,6,7,4,2)

find_peaks(x = test_1, ignore_threshold = 0, span = 3, strict = TRUE)



