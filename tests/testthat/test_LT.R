rm(list = ls())

F_mx <- HMD3mx$female$SWE
ages <- as.numeric(rownames(F_mx))
year <- 1970
mx <- F_mx[, paste(year)]

LT <- lifetable(x = ages, mx = mx, sex = "female")
LT
LT$lt

test_that("Life table tests", {
  expect_identical(class(LT), "lifetable")
  expect_identical(class(LT$lt), "data.frame")
  expect_false(is.null(LT))
  expect_true(all(LT$lt >= 0))
})


# Example 2 --- Abridge life table ------------
x  = c(0, 1, seq(5, 110, by = 5))
mx = c(.053, .005, .001, .0012, .0018, .002, .003, .004, 
       .004, .005, .006, .0093, .0129, .019, .031, .049, 
       .084, .129, .180, .2354, .3085, .390, .478, .551)
lt = lifetable(x, mx = mx, sex = "female")
lt
