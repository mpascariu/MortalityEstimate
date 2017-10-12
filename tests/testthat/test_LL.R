rm(list = ls())
library(MortalityEstimate)
library(testthat)

F_mx <- HMD3mx$female$SWE
ages <- as.numeric(rownames(F_mx))
year <- 1970
mx <- F_mx[, paste(year)]

LT <- lifetable(x = ages, mx = mx, sex = "female")
LT


test_that("Life table tests", {
  expect_identical(class(LT), "lifetable")
  expect_identical(class(LT$lt), "data.frame")
  expect_true(all(LT$lt >= 0))
})



# Select the 1965 - 1990 time interval and fit the Linear-Link model
ages  <- 0:100 # available ages in our datasets
years <- 1965:1990 # available years
sex   <- 'female'
SWEmx <- HMD3mx$female$SWE[paste(ages), paste(years)]

fit_LL <- LinearLink(mx = SWEmx,
                     mx_ages = ages,
                     mx_years = years,
                     mx_country = 'SWEDEN',
                     theta = 0,
                     method = 'LSE')
fit_LL

test_that("LinearLink Tests", {
  expect_true(!is.null(fit_LL))
})


