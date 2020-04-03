% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCVAR_estn.R
\name{FCVARestn}
\alias{FCVARestn}
\title{Estimate FCVAR model}
\usage{
FCVARestn(x, k, r, opt)
}
\arguments{
\item{x}{A matrix of variables to be included in the system.}

\item{k}{The number of lags in the system.}

\item{r}{The cointegrating rank.}

\item{opt}{A list object that stores the chosen estimation options,
generated from \code{EstOptions()}.}
}
\value{
A list object \code{results} containing the estimation results,
including the following parameters:
\describe{
  \item{\code{startVals}}{Starting values used for optimization.}
  \item{\code{options}}{Estimation options.}
  \item{\code{like}}{Model log-likelihood.}
  \item{\code{coeffs}}{Parameter estimates.}
  \item{\code{rankJ}}{Rank of Jacobian for the identification condition.}
  \item{\code{fp}}{Number of free parameters.}
  \item{\code{SE}}{Standard errors.}
  \item{\code{NegInvHessian}}{Negative of inverse Hessian matrix.}
  \item{\code{Residuals}}{Model residuals.}
  \item{\code{cPolyRoots}}{Roots of characteristic polynomial.}
}
}
\description{
\code{FCVARestn} estimates the Fractionally Cointegrated VAR model.
  It is the central function in the \code{FCVAR} package with several nested functions, each
	described below. It estimates the model parameters, calculates the
	standard errors and the number of free parameters, obtains the residuals
	and the roots of the characteristic polynomial, and prints the output.
}
\examples{
opt <- EstOptions()
x <- data(JNP2014)
FCVARestn(x,k = 3,r = 1,opt)
}
\seealso{
\code{EstOptions} to set default estimation options.
\code{FCVARestn} calls this function at the start of each estimation to verify
validity of options.
}
\concept{FCVAR estimation functions}