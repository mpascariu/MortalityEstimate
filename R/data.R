# Data in the package

#' Death rates for female and male populations in US, France and Sweden 
#'
#' Dataset containing list containing 6 matrices with death rates for 
#' female and male populations in  USA, France and Sweden between 1965 and 2014. 
#' This data is provided for testing purposes only. It may not be up to date anymore.
#' Download the actual data free of charge from \url{http://www.mortality.org}.  
#'
#' @format A list with 6 matrices with 101 rows and 50 columns:
#' \describe{
#'   \item{rows}{age interval}
#'   \item{colums}{years}
#' }
#' @source Human Mortality Database, \url{http://www.mortality.org}.
"HMD.test.data"


#' 719 life tables from HMD (2012). 
#' 
#' This data is provided for testing purposes only. 
#' This data is provided for testing purposes only. It may not be up to date anymore.
#' Download the actual data free of charge from \url{http://www.mortality.org}.
"HMD719"


#' Import Packages and functions
#'
#' @importFrom stats complete.cases lsfit optim 
#' smooth.spline coef reshape fitted poisson median predict uniroot
#' @importFrom pbapply startpb closepb setpb
NULL