# Data in the package

#' Death rates for female and male populations in US, France and Sweden 
#'
#' Dataset containing list containing 6 matrices with death rates for 
#' female and male populations in USA, France and Sweden between 1965 and 2014. 
#' This data is provided for testing purposes only. It may not be up to date anymore.
#' Download the actual data free of charge from \url{http://www.mortality.org}.  
#' @source Human Mortality Database, \url{http://www.mortality.org}.
"HMD3mx"


#' 719 life tables from HMD 
#' 
#' Data used in the Wilmoth et. al. (2012) article. Today this dataset is outdated,
#' download the actual data free of charge from \url{http://www.mortality.org}. 
#' @references John Wilmoth, Sarah Zureick, Vladimir Canudas-Romo, Mie Inoue & 
#' Cheryl Sawyer (2012): A flexible two-dimensional mortality model for use in 
#' indirect estimation, Population Studies: A Journal of Demography, 66:1, 1-28
#' \url{http://dx.doi.org/10.1080/00324728.2011.611411}
#' @source Human Mortality Database, \url{http://www.mortality.org}. (2012)
"HMD719"


#' Import Packages and functions
#'
#' @importFrom stats complete.cases lsfit optim 
#' smooth.spline coef reshape fitted poisson median predict uniroot
#' @importFrom pbapply startpb closepb setpb
#' @name foo_imports
#' @keywords internal
NULL
