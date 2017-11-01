# rm(list = ls())

# DATA
sex = "female"
HMD719f <- HMD719[HMD719$sex == sex, ]

# Fit model
x <- c(0,1, seq(5, 110, by = 5))
W <- wilmoth(x, LT = HMD719f, sex = sex)
W

test.wilmoth <- function(Y) {
  expect_true(class(Y) == "wilmoth")
  expect_false(is.null(Y))
  expect_false(is.null(coef(Y)))
  expect_false(any(is.na(fitted(Y))))
  expect_false(any(is.na(resid(Y))))
  expect_false(any(is.na(coef(Y))))
}

test.wilmoth(W)

# Build life tables with various choices of 2 input parameters ---
L1 <- wilmothLT(W, q0_5 = 0.05, k = 0.1)
L2 <- wilmothLT(W, q0_5 = 0.05, e0 = 65)
L3 <- wilmothLT(W, q0_5 = 0.05, q15_45 = 0.2)
L4 <- wilmothLT(W, q0_5 = 0.05, q15_35 = 0.125)
L5 <- wilmothLT(W, q0_1 = 0.01, k = 0.1)
L6 <- wilmothLT(W, q0_1 = 0.01, e0 = 65)
L7 <- wilmothLT(W, q0_1 = 0.05, q15_45 = 0.2)
L8 <- wilmothLT(W, q0_1 = 0.05, q15_35 = 0.125)
L9 <- wilmothLT(W, k = 0.01, e0 = 65)
L10 <- wilmothLT(W, k = 0.01, q15_45 = 0.2)
L11 <- wilmothLT(W, k = 0.01, q15_35 = 0.125)
L12 <- wilmothLT(W, q15_45 = 0.125, e0 = 65)


test.wilmothLT <- function(P){
  expect_true(class(P) == "wilmothLT")
  expect_false(is.null(P))
  expect_false(any(is.na(P$lt)))
  expect_true(all(P$lt[,-1] >= 0))
  expect_false(any(is.na(P$values)))
  
}

for (i in 1:12) test.wilmothLT(get(paste0("L",i)))

# Test Q5 ---
Q5 = 1 - (1 - L1$lt$qx[1])*(1 - L1$lt$qx[2])
expect_equal(Q5, 0.05)
