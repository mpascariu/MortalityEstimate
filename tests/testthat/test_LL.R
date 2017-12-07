rm(list = ls())

ages  <- 0:100 
years <- 1965:1990
sex   <- 'female'
SWEmx <- HMD4mx$SWE[paste(ages), paste(years)]

M1 <- LinearLink(x = ages, mx = SWEmx, y  = years,
                 country = 'SWEDEN', theta   = 0, method = 'LSE')

M2 <- LinearLink(x = ages, mx = SWEmx, y = years,
                 country = 'SWEDEN', theta = 0, method = 'MLE')


test.LinearLink <- function(X){
  expect_false(is.null(X))
  expect_false(any(is.na(fitted(X))))
  expect_false(any(is.na(resid(X))))
  expect_false(any(is.na(X$fitted.life.tables)))
  expect_false(any(is.na(coef(X))))
  expect_false(is.null(summary(X)))
  expect_false(is.null(coef(X)))
}

test.LinearLink(M1)
test.LinearLink(M2)

M1$fitted.life.tables[is.na(M1$fitted.life.tables$ex), ]

M1$fitted.life.tables[99:101, ]


P1 <- LinearLinkLT(M1, ex = 85)
P2 <- LinearLinkLT(M2, ex = 85)


test.LinearLinkLT <- function(P){
  expect_false(is.null(P))
  expect_false(any(is.na(P)))
  expect_false(any(is.na(P$lt)))
  expect_true(all(P$lt[,-1] >= 0))
}

test.LinearLinkLT(P1)
test.LinearLinkLT(P2)
