#install_github("mpascariu/LinearLink")
#library(mnormt)
#library(devtools)
library(latticeExtra)
library(dplyr)
library(LinearLink)

path <- paste0(getwd(), '/data-raw/LinearLInk MLE/')

Deaths    <- as.matrix(read.table(paste0(path,"SWIdeaths.txt"),T))[,16:41]
Exposures <- as.matrix(read.table(paste0(path,"SWIexposures.txt"),T))[,16:41]
mx <- Deaths/Exposures

# # Example from Pascariu et al 2016 --------------------------------------

# Select the 1965 - 1990 time interval and fit the Linear-Link model
ages    <- 0:100 # available ages in our datasets
years   <- 1965:1990 # available years
fit_SWE <- LinearLink(mx = mx, ages, years, 'SWEDEN', x_fit = 0, method = 'LSE')
summary(fit_SWE)

# Derive a mortality curve (life table) from a value of
# life expectancy at birth in 2014, say 84.05
new_e0   <- 84.05
pred_SWE <- predict(fit_SWE, new_e0)
pred_SWE$lt

Coef1 <- list(fit_SWE$coefficients,fit_SWE$k_values)

# Fitting model with ML ---------------------------------------------------

# Based mainly in Brouhns et al 2002
# Some useful functions to fit the model 
# We just reparametrize to fit Pascarius et al's model
# Alpha= b_x log(e(theta))
# Beta= v_x
# Kappa= k
# Note we need Deaths an Exposures

## Update Alpha

Update.alpha <- function(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit){
  difD <- Dth - D.fit
  Alpha <- Alpha + difD %*% One / (D.fit %*% One)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Alpha = Alpha, D.fit = D.fit)
}

## Update Beta
Update.beta <- function(Alpha, Beta, Kappa,
                        One, Dth, Exp, D.fit){
  difD <- Dth - D.fit
  Kappa2 <- Kappa * Kappa
  Beta <- Beta + difD %*% Kappa / (D.fit %*% Kappa2)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Beta = Beta, D.fit = D.fit)
}

## Update Kappa
Update.kappa <- function(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit){
  difD <- Dth - D.fit
  Beta2 <- Beta * Beta
  Kappa <- Kappa + t(difD) %*% Beta / (t(D.fit) %*% Beta2)
  Kappa <- Kappa - mean(Kappa)
  Kappa <- Kappa / sqrt(sum(Kappa * Kappa))
  Kappa <- matrix(Kappa, ncol = 1)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Kappa = Kappa, D.fit = D.fit)
}

##
PoissonMl_Estima <- function(Dth, Exp){
  # dimensions
  n <- ncol(Dth)
  # Initialise
  One <- matrix(1, nrow = n, ncol = 1)    
  Fit.init <- log((Dth + 1)/(Exp + 2))
  Alpha <- Fit.init %*% One / n
  Beta <- matrix(1 * Alpha, ncol = 1)
  sum.Beta <- sum(Beta) 
  Beta <- Beta / sum.Beta
  Kappa <- matrix(seq(n, 1, by = -1), nrow = n, ncol = 1)
  Kappa <- Kappa - mean(Kappa)
  Kappa <- Kappa / sqrt(sum(Kappa * Kappa))
  Kappa <- Kappa * sum.Beta
  
  # Iteration
  D.fit <- Exp * exp(Fit.init)
  for (iter in 1:50) {
    Alpha.old <- Alpha
    Beta.old <- Beta
    Kappa.old <- Kappa
    #
    temp <- Update.alpha(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit)
    D.fit <- temp$D.fit
    Alpha <- temp$Alpha
    #
    temp <- Update.beta(Alpha, Beta, Kappa,
                        One, Dth, Exp, D.fit)
    D.fit <- temp$D.fit
    Beta <- temp$Beta
    #
    temp <- Update.kappa(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit)
    D.fit <- temp$D.fit
    Kappa <- temp$Kappa
    crit <- max(max(abs(Alpha - Alpha.old)),
                max(abs(Beta - Beta.old)),
                max(abs(Kappa - Kappa.old)))
    if (crit < 1e-04) break
  }
  # constraints
  sum.Beta <- sum(Beta)
  Beta <- Beta / sum.Beta
  Kappa <- Kappa * sum.Beta
  # output
  out <- list(Alpha = Alpha, Beta = Beta, Kappa = Kappa)
  return(out)
}



Model.ML <- PoissonMl_Estima(Dth = Deaths, Exp = Exposures)

# Plot parameters ---------------------------------------------------

a.hat       <- Model.ML$Alpha
b.hat       <- Model.ML$Beta
k.hat       <- Model.ML$Kappa

windows(record = T)
par(mfrow = c(1,3))
plot(coef(fit_SWE)$bx, ylim = c(-8.5, 0))
lines(a.hat, lwd = 2, col = 2)

plot(coef(fit_SWE)$vx, ylim = c(-0.1, 0.3))
lines(b.hat, lwd = 2, col = 2)

plot(fit_SWE$k_values$k, ylim = c(-25, 25))
lines(k.hat, lwd = 2, col = 2)

One <- matrix(1, nrow = length(years), ncol = 1)

# estimated rates
MU.hat      <- exp(a.hat %*% t(One) + b.hat %*% t(k.hat))
# or in your case k=0 to optimize

plot(log(fit_SWE$fitted.values[, 10]))
lines(log(MU.hat[, 10]))

# then just follow the smoothing and the optimization of the curve

#Comparison with Pascarius' estimation
f1 <- xyplot(a.hat ~ ages, col = "black", xlab = "Age", 
             main = "a_x = b_xlog(e(theta))",type="l",
             panel = function(x, y, ...){   
               panel.grid(h = -1,v = 0,col = 'dark grey', lty = 3)
               panel.abline(v = c(seq(0,100,10)),col = 'dark grey', lty = 3)       
               panel.xyplot(x, y, ...)
             }
)
f2 <- xyplot(b.hat ~ ages, col = "black", xlab = "Age", main = "v_x", type = "l",
             panel = function(x, y, ...){   
               panel.grid(h = -1, v = 0, col = 'dark grey', lty = 3)
               panel.abline(v = c(seq(0,100,10)), col = 'dark grey', lty = 3)       
               panel.xyplot(x, y, ...)
             }
)
f3  <- xyplot(k.hat ~ 1965:1990, col="black",xlab="Year", main="k",type="l",
              panel = function(x, y, ...){   
                panel.grid(h=-1,v=0,col='dark grey',lty=3)
                panel.abline(v=c(seq(1950,2007,10)),col='dark grey',lty=3)       
                panel.xyplot(x, y, ...)
              }
)


require(gridExtra)
grid.arrange(f1, f2,f3, nrow=1)
dev.off()

# --------------------------------------

library(StMoMo)
lc

#sum(kt) = 0 and log link
LC1 <- lc()
LCfit1 <- fit(LC1, Dxt = EWMaleData$Dxt,Ext = EWMaleData$Ext, 
            ages = EWMaleData$ages, years = EWMaleData$years,
            ages.fit = 55:89)
plot(LCfit1)



