# rm(list = ls())

ages  <- 0:100 
years <- 1965:1990
sex   <- 'female'
SWEmx <- HMD3mx$female$SWE[paste(ages), paste(years)]

mdl <- LinearLink(x = ages,
                  mx = SWEmx,
                  y = years,
                  country = 'SWEDEN',
                  theta = 0,
                  method = 'LSE')
mdl

test_that("LinearLink Tests", {
  expect_false(is.null(mdl))
  expect_false(is.null(summary(mdl)))
  expect_false(is.null(coef(mdl)))
})


W2 = wilmoth(x = ages, mx = SWEmx, sex = sex)

