# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# Last update: Fri Apr 30 11:41:40 2021
# --------------------------------------------------- #
remove(list = ls())

ages  <- 0:100 
years <- 1965:1990
sex   <- 'female'
SWEmx <- HMD4mx$SWE[paste(ages), paste(years)]

M1 <- LinearLink(x = ages, 
                 mx = SWEmx, 
                 y  = years,
                 country = 'SWEDEN', 
                 theta   = 0, 
                 method = 'LSE')

M2 <- LinearLink(x = ages, 
                 mx = SWEmx, 
                 y = years,
                 country = 'SWEDEN', 
                 theta = 0, 
                 method = 'MLE')

M3 <- LinearLink(x = ages, 
                 mx = SWEmx, 
                 y = years,
                 country = 'SWEDEN', 
                 theta = 10, 
                 method = 'MLE')


test.LinearLink <- function(X){
  expect_false(is.null(X))
  expect_false(any(is.na(fitted(X))))
  expect_false(any(is.na(resid(X))))
  expect_false(any(is.na(X$fitted.life.tables)))
  expect_false(any(is.na(coef(X))))
  expect_false(is.null(summary(X)))
  expect_false(is.null(coef(X)))
  expect_output(print(X))
  expect_output(print(summary(X)))
}

test.LinearLink(M1)
test.LinearLink(M2)

M1$fitted.life.tables[is.na(M1$fitted.life.tables$ex), ]

M1$fitted.life.tables[99:101, ]


P1 <- LinearLinkLT(M1, ex = 85)
P2 <- LinearLinkLT(M2, ex = 85)
P3 <- LinearLinkLT(M1, ex = 85, use.vx.rotation = TRUE)
P4 <- LinearLinkLT(M1, ex = 85, use.vx.rotation = TRUE, e0_threshold = 90)
P5 <- LinearLinkLT(M1, ex = 85, use.vx.rotation = TRUE, e0_u = 80)

test.LinearLinkLT <- function(P){
  expect_false(is.null(P))
  expect_false(any(is.na(P)))
  expect_false(any(is.na(P$lt)))
  expect_true(all(P$lt[,-1] >= 0))
}

test.LinearLinkLT(P1)
test.LinearLinkLT(P2)
test.LinearLinkLT(P3)
test.LinearLinkLT(P4)
test.LinearLinkLT(P5)

# ----------------------------------------------
# Test Data
expect_false(is.null(HMD4mx))
expect_output(print.MortalityEstimateData(HMD4mx))


# ----------------------------------------------------------------------------
# Test messages

expect_message(
  LinearLink(x = ages,
             mx = SWEmx,
             y  = years,
             country = 'SWEDEN',
             theta   = 51,
             method = 'LSE')
)

expect_error(
  LinearLink(x = ages[-1],
             mx = SWEmx,
             y  = years,
             country = 'SWEDEN',
             theta   = 0,
             method = 'LSE')
)

expect_error(
  LinearLink(x = ages,
             mx = SWEmx,
             y  = years[-1],
             country = 'SWEDEN',
             theta   = 0,
             method = 'LSE')
)

expect_error(
  LinearLink(x = ages,
             mx = SWEmx,
             y  = years,
             country = 'SWEDEN',
             theta   = 0,
             method = 'some_typo')
)

expect_error(
  LinearLinkLT(M3, ex = 85, use.vx.rotation = TRUE)
)









