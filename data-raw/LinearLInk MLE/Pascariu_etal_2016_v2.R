#install_github("mpascariu/LinearLink")
#library(mnormt)
#library(devtools)
# setwd("C:/Users/jmaburto/Desktop/Pascariu_etal2016LL/ML model")
rm(list = ls())
library(latticeExtra)
library(dplyr)
library(LinearLink)
library(reshape2)

Deaths      <- as.matrix(read.table("SWIdeaths.txt",T))[,16:41]
Exposures   <- as.matrix(read.table("SWIexposures.txt",T))[,16:41]
mx <- Deaths/Exposures

# # Example from Pascariu et al 2016 --------------------------------------

# Select the 1965 - 1990 time interval and fit the Linear-Link model
ages    <- 0:100 # available ages in our datasets
years   <- 1965:1990 # available years
fit_SWE <- LinearLink(mx = mx, mx_ages = ages, mx_years = years,
                      mx_country = 'SWEDEN', theta = 0, method = 'LSE',use.smooth = F)
plot(fit_SWE$coefficients$bx)

# Derive a mortality curve (life table) from a value of
# life expectancy at birth in 2014, say 84.05
new_e0   <- 84.05
pred_SWE <- predict(fit_SWE, new_e0)
pred_SWE$lt

Coef1 <- list(fit_SWE$coefficients,fit_SWE$k_values)
sum(Coef1[[1]]$bx)
sum(Coef1[[1]]$vx)
sum(fit_SWE$k_values$k)

# Fitting model with ML ---------------------------------------------------

#matrix with life expectancies
ex_theta <- acast(fit_SWE$fitted.life.tables, age~year, value.var = 'ex')
Model.ML <- PMLE(Dx = Deaths, Ex = Exposures, log_ex_theta = log(ex_theta[1,]))


# Comparison
bx1 <- coef(fit_SWE)$bx
vx1 <- coef(fit_SWE)$vx
k1 <- coef(fit_SWE)$k

bx2 <- Model.ML$bx
vx2 <- Model.ML$vx
k2 <- Model.ML$k
plot(vx2)

par(mfrow = c(1, 4))
plot(bx1, col = 4, pch = 16, main = 'bx')
lines(bx2, col = 2, lwd = 2)

plot(vx1, col = 4, pch = 16, main = 'vx')
lines(vx2, col = 2, lwd = 2)

plot(k1, col = 4, pch = 16, main = 'k')
lines(k2, col = 2, lwd = 2)

plot(log(mx[,1]), pch = 16, main = 'mx')
lines(log(fit_SWE$fitted.values[, '1965']), col = 4, lwd = 2)
lines(bx2*log(75.2) + k2[1]*vx2, col = 2, lwd = 2)





