#' A package for estimating the Fractionally Cointegrated VAR model.
#'
#' The FCVAR package estimates the Fractionally Cointegrated
#' Vector Autoregressive (VAR) model. It includes functions for lag selection,
#' cointegration rank selection and hypothesis testing.
#'
#' Functions in the FCVAR package are divided into four categories:
#' Estimation, Postestimation, Specification and Auxiliary functions.
#'
#' @section Estimation functions:
#' The estimation functions include the primary estimation function \code{FCVARestn}
#' and associated functions to set estimation options and display results.
#' Some of these functions define, modify and test the user-specified options for estimation.
#' \code{FCVARoptions} defines the default estimation options used in the FCVAR
#' estimation procedure and the related programs.
#' The user can then revise the options such as the settings for optimization and
#' restrictions for testing hypotheses.
#' After making these changes, an internal function \code{FCVARoptionUpdates} sets and tests
#' estimation options for validity and compatibility.
#'
#' @section Postestimation functions:
#' The postestimation functions are used to display summary statistics,
#' test hypotheses and test the goodness of fit of the estimated model.
#'   These include:
#' \describe{
#'   \item{\code{FCVARhypoTest}}{for a likelihood ratio test of a restricted
#'   vs. an unrestricted model}
#'   \item{\code{FCVARboot}}{for generating a distribution of a likelihood ratio
#'   test statistic}
#'   \item{\code{FCVARforecast}}{for calculating recursive forecasts with the FCVAR model}
#' }
#'
#' @section Specification functions:
#' The specification functions are used to estimate a series of models in order
#'   to make model specfication decisions.
#'   These include:
#' \describe{
#'   \item{\code{FCVARlagSelect}}{for selection of the lag order}
#'   \item{\code{FCVARrankTests}}{for choosing the cointegrating rank}
#'   \item{\code{FCVARbootRank}}{for generating a distribution of a likelihood ratio
#'  test statistic for the rank test}
#' }
#'
#' @section Auxiliary functions:
#' The auxiliary functions are used to perform intermediate calculations
#'   for estimation.
#'   These functions are mainly designed for use only within the estimation function.
#'   Some exceptions include:
#' \describe{
#'   \item{\code{FracDiff}}{for fractionally differencing a multivariate series}
#'   \item{\code{FCVARsimBS}}{for generating bootstrap samples from the FCVAR model}
#'   \item{\code{FCVARlikeGrid}}{for performing a grid-search optimization with the FCVAR likelihood function}
#' }
#'
#' @section Examples:
#' A dataset \code{votingJNP2014} is included for examples of the model building process.
#'   Sample model builds with hypothesis tests and examples of other extensions are found
#'   in the example script \code{FCVAR_demo_JNP2014.R}.
#'   See FCVAR_README.pdf for details
#'   at
#'
#'   \code{https://github.com/LeeMorinUCF/FCVAR/blob/master/FCVAR_README.pdf}
#'
#'   and also see \code{https://sites.google.com/view/mortennielsen/software}
#'   for more information about estimating the FCVAR model.
#'
#' @docType package
#' @name FCVAR
#' @return Returns \code{NULL}. Object included for description only.
NULL


#' Aggregate support for Canadian political parties.
#'
#' A dataset containing the aggregate support for Canadian political
#' parties and economic indicators from Canada and the United States.
#'
#' @format A data frame with 316 rows and 6 variables:
#' \describe{
#'   \item{lib}{aggregate support for the Liberal party}
#'   \item{pc}{aggregate support for the Conservative party}
#'   \item{ir_can}{Canadian 3-month T-bill rates}
#'   \item{ir_us}{US 3-month T-bill rates}
#'   \item{un_can}{Canadian unemployment rate}
#'   \item{un_us}{US unemployment rate}
#' }
#' @source \url{https://sites.google.com/view/mortennielsen/software}
"votingJNP2014"

