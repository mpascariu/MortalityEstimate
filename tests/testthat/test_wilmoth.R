rm(list = ls())

# DATA
sex = "female"
HMD719f <- HMD719[HMD719$sex == sex, ]

# Fit model
W <- wilmoth(data = HMD719f, sex)

W

# Build life tables with various choices of 2 input parameters ---
# (case 1) Using 5q0 and k
WLT_1 <- wilmothLT(W, q0_5 = 0.05, k = 0.1)
  
test_that("wilmoth tests", {
  expect_false(is.null(W))
  expect_false(is.null(coef(W)))
})