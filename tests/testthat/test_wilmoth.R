# rm(list = ls())

# DATA
sex = "female"
HMD719f <- HMD719[HMD719$sex == sex, ]

# Fit model
x <- c(0,1, seq(5, 110, by = 5))
W <- wilmoth(x, LT = HMD719f, sex = sex)



# Build life tables with various choices of 2 input parameters ---
# (case 1) Using 5q0 and k
L1 <- wilmothLT(W, q0_5 = 0.05, k = 0.1)

Q5 = 1 - (1 - L1$lt$qx[1])*(1 - L1$lt$qx[2])

test_that("wilmoth tests", {
  expect_false(is.null(W))
  expect_false(is.null(coef(W)))
  expect_equal(Q5, 0.05)
})

