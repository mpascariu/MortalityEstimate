rm(list = ls())
# Install latest version of LinearLink package 
# devtools::install_github('mpascariu/LinearLink')
library(LinearLink)
library(gnm)

# Select the 1965 - 1990 time interval and fit the Linear-Link model
ages  <- 0:100 # available ages in our datasets
years <- 1965:1990 # available years

SWEmx <- HMD.test.data$SWE[paste(ages), paste(years)]

fit_LL <- LinearLink(mx = SWEmx, mx_ages = ages, mx_years = years, 
                     mx_country = 'SWEDEN', theta = 0, method = 'LSE')
fit_LL

summary(fit_LL) # summary
coef(fit_LL) # cofficients

bx1 <- coef(fit_LL)$bx
vx1 <- coef(fit_LL)$vx
k1 <- coef(fit_LL)$k

# ---------------------------------
# Create a data.frame for the MLE model
dta <- fit_LL$fitted.life.tables 
dta$ex0 <- NA
dta$Ex <- 1

# create a colums with life expectancy at birth in every year (over all ages)
yr <- unique(dta$year)
for (i in 1:yr) {
  ex0_i <- dta[dta$year == yr[i] & dta$age == 0, 'ex']
  dta[dta$year == yr[i], ]$ex0 <- ex0_i
}
# I'm not sure why do I get errors. But it works.
head(dta)
tail(dta)

# w poisson MLE
set.seed(8)
LL <- gnm(dx ~ -1 + as.factor(age):log(ex0) +
          Mult(as.factor(age), as.factor(year)) +
          offset(log(lx)), family = poisson(link = "log"), data = dta)

# w binomial MLE
# set.seed(10)
# LL <- gnm(dx/lx ~ -1 + as.factor(age):log(ex0) + Mult(as.factor(age), as.factor(year)),
#                     weights = lx, family = binomial(link = "logit"), data = dta)

# w normal MLE
# LL <- gnm(dx/lx ~ -1 + as.factor(age):log(ex0) + Mult(as.factor(age), as.factor(year)),
                # family = gaussian(link="log"), data = dta)

LLcoef <- data.frame(coef(LL))
bx2 <- LLcoef[1:101, ]
vx2 <- LLcoef[102:202, ]
k2 <- LLcoef[203:228, ]


mx1965 <- log(HMD.test.data$SWE[, '1965'])
mx1965.LL1 <- log(fit_LL$fitted.values[, '1965'])
mx1965.LL2 <- bx2*log(dta[dta$year == 1965,'ex0'][1]) + k2[1]*vx2

# Comparison
par(mfrow = c(1, 4))
plot(bx1, col = 4, pch = 16, main = 'bx')
lines(bx2, col = 2, lwd = 2)
plot(vx1, col = 4, pch = 16, main = 'vx')
lines(vx2, col = 2, lwd = 2)

plot(k1, col = 4, pch = 16, main = 'k')
lines(k2, col = 2, lwd = 2)
plot(mx1965, pch = 16, main = 'mx')
lines(mx1965.LL1, col = 4, lwd = 2)
lines(mx1965.LL2, col = 2, lwd = 2)






