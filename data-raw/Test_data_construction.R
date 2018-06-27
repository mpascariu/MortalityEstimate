rm(list = ls())
library(dplyr)
library(MortalityLaws)

#######################  Setup  #################################
# Test data to be added in the package
# Import HMD dataset (used for fitting the model)

cntr = c('SWE', 'FRATNP', 'USA')
HMDuser = "username@email.com"
HMDpass = "HMDpassword"

HMD_Dx <- ReadHMD(what = "Dx", countries = cntr,
                  username = HMDuser, password = HMDpass,
                  save = F)$data
HMD_Ex <- ReadHMD(what = "Ex", countries = cntr,
                  username = HMDuser, password = HMDpass,
                  save = F)$data



build.mx.wide <- function(HMD.Ex, HMD.Dx, years, ages, sex) {
  Ext <- HMD.Ex %>% select(country:Age, get(sex)) %>% filter(Year %in% years)
  Dxt <- HMD.Dx %>% select(country:Age, get(sex)) %>% filter(Year %in% years)
  colnames(Ext)[4] <- "Ex"
  colnames(Dxt)[4] <- "Dx"
  mxt  <- left_join(Ext, Dxt, by = c('country', 'Year','Age')) %>% mutate(mx = Dx/Ex)
  cntr <- levels(HMD_Dx$country)
  
  HMD_mx <- list()
  for (i in 1:length(cntr)) {
       cntr_i <- cntr[i]
       mx_i <- mxt[complete.cases(mxt), ] %>% 
              filter(Year %in% year.range, country == cntr_i, Age %in% ages) %>% 
              select(country:Age, mx) %>% 
              reshape(., direction = 'wide', idvar = c('country','Age'), timevar = 'Year') %>% 
              select(-(country:Age))
       dimnames(mx_i) <- list(ages, years)
       HMD_mx[[i]] <- mx_i
  }
  names(HMD_mx) <- cntr
  return(HMD_mx)
}



age.range  <- 0:100
year.range <- 1965:2014
HMD_mxm <- build.mx.wide(HMD_Ex, HMD_Dx, year.range, age.range, sex = "Male")
HMD_mxf <- build.mx.wide(HMD_Ex, HMD_Dx, year.range, age.range, sex = "Female")

HMD3mx = list(female = HMD_mxf, male = HMD_mxm)
devtools::use_data(HMD3mx, overwrite = TRUE)
#----------------------------------------------------------------------------






