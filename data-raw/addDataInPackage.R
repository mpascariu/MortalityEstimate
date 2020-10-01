# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: MIT
# Last update: Wed Sep 30 10:02:26 2020
# --------------------------------------------------- #
remove(list = ls())
library(dplyr)
library(MortalityLaws)

#######################  Setup  #################################
# Test data to be added in the package
# Import HMD dataset (used for fitting the model)

cntr = c("GBRTENW", "FRATNP", "SWE", "USA")
HMDuser = "..."
HMDpass = "..."

HMD_Dx <- ReadHMD(
            what = "Dx", 
            countries = cntr,
            username = HMDuser, 
            password = HMDpass,
            save = FALSE)$data[, 1:4]

HMD_Ex <- ReadHMD(
            what = "Ex", 
            countries = cntr,
            username = HMDuser, 
            password = HMDpass,
            save = F)$data[, 1:4]

HMD_mx <- data.frame(HMD_Dx[, 1:3], mx = HMD_Dx[, 4]/HMD_Ex[, 4])

build.mx.wide <- function(HMD_mx, years, ages) {
  HMDmx <- HMD_mx %>% filter(Year %in% years, Age %in% ages)
  out   <- list()
  for (i in 1:length(cntr)) {
    mx_i <- HMDmx %>% 
      filter(country == cntr[i]) %>% 
      reshape(
        direction = 'wide', 
        idvar = c('country','Age'), 
        timevar = 'Year') %>% 
      select(-(country:Age))
    dimnames(mx_i) <- list(ages, years)
    out[[i]] <- mx_i
  }
  names(out) <- cntr
  return(out)
}


years = 1965:2018
ages  = 0:110
dta  <- build.mx.wide(HMD_mx, years, ages)
HMD4mx <- structure(class = "MortalityEstimateData", dta)

usethis::use_data(HMD4mx, overwrite = TRUE)
#----------------------------------------------------------------------------






