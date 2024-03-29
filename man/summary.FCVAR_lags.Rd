% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCVAR_spec.R
\name{summary.FCVAR_lags}
\alias{summary.FCVAR_lags}
\title{Summarize Statistics from Lag Order Selection}
\usage{
\method{summary}{FCVAR_lags}(object, ...)
}
\arguments{
\item{object}{An S3 object of type \code{FCVAR_lags} containing the results
from repeated estimation of the FCVAR model with different orders
of the autoregressive lag length. It is the output of \code{FCVARlagSelect}.}

\item{...}{additional arguments affecting the summary produced.}
}
\description{
\code{summary.FCVAR_lags} prints a summary of the table of statistics from
the output of \code{FCVARlagSelect}.
\code{FCVARlagSelect} takes a matrix of variables and performs lag
	selection on it by using the likelihood ratio test.
}
\examples{
\donttest{
opt <- FCVARoptions()
opt$gridSearch   <- 0 # Disable grid search in optimization.
opt$dbMin        <- c(0.01, 0.01) # Set lower bound for d,b.
opt$dbMax        <- c(2.00, 2.00) # Set upper bound for d,b.
opt$constrained  <- 0 # Impose restriction dbMax >= d >= b >= dbMin ? 1 <- yes, 0 <- no.
x <- votingJNP2014[, c("lib", "ir_can", "un_can")]
FCVAR_lag_1 <- FCVARlagSelect(x, kmax = 3, r = 3, order = 12, opt)
summary(object = FCVAR_lag_1)
}
}
\seealso{
\code{FCVARoptions} to set default estimation options.
\code{FCVARestn} is called repeatedly within this function
for each candidate lag order.
\code{summary.FCVAR_lags} prints a summary of the output of \code{FCVARlagSelect} to screen.

Other FCVAR specification functions: 
\code{\link{FCVARbootRank}()},
\code{\link{FCVARlagSelect}()},
\code{\link{FCVARrankTests}()},
\code{\link{summary.FCVAR_ranks}()}
}
\concept{FCVAR specification functions}
