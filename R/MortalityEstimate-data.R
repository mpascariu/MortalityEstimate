# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Sat Dec  1 21:20:51 2018
# --------------------------------------------------- #

# Data in the package ----


#' Death Rates for Female Populations in England and Wales, France, Sweden and USA 
#'
#' Dataset containing containing 4 matrices with death rates for 
#' female populations in England and Wales, France, Sweden and USA 
#' between 1965 and 2014. This data is provided for testing purposes only. 
#' It may not be up to date anymore. Download the actual data free of charge 
#' from \url{http://www.mortality.org}.  
#' @source Human Mortality Database, \url{http://www.mortality.org}.
"HMD4mx"

#' Print function for HMD4mx data
#' 
#' @param x A \code{MortalityEstimateData} object.
#' @param ... Further arguments passed to or from other methods.
#' @export
#' @keywords internal
print.MortalityEstimateData <- function(x, ...) {
  cat("\nMortalityEstimate Test Data\n")
  cat(" Populations : England and Wales, France, Sweden, USA\n")
  cat(" Series      : Death rates for females (raw)\n")
  cat(" Years       : 1965 - 2014\n")
  cat(" Ages        : 0 - 110\n")
  cat(" Format      : List containg 4 data frames\n")
  cat(" Source      : Human Mortality Database\n")
  cat(" Download    : December 7, 2017")
}


#' 719 Life Tables from HMD 
#' 
#' Data used in the Wilmoth et. al. (2012) article. Today this dataset is outdated,
#' download the actual data free of charge from \url{http://www.mortality.org}. 
#' @references John Wilmoth, Sarah Zureick, Vladimir Canudas-Romo, Mie Inoue & 
#' Cheryl Sawyer (2012): A flexible two-dimensional mortality model for use in 
#' indirect estimation, Population Studies: A Journal of Demography, 66:1, 1-28
#' \url{http://dx.doi.org/10.1080/00324728.2011.611411}
#' @source Human Mortality Database, \url{http://www.mortality.org}. (2012)
"HMD719"

