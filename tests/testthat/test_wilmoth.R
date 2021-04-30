# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# Last update: Fri Apr 30 13:05:38 2021
# --------------------------------------------------- #
remove(list = ls())


# DATA
x <- c(0,1, seq(5, 110, by = 5))
H0 <- HMD719[HMD719$sex == "female", ]
H1 <- array(H0$mx, dim = c(length(x), nrow(H0)/length(x)))

# Fit model
x <- c(0,1, seq(5, 110, by = 5))
W0 <- wilmoth(x, LT = H0)
W1 <- wilmoth(x, mx = H1)


test.wilmoth <- function(Y) {
  expect_true(class(Y) == "wilmoth")
  expect_false(is.null(Y))
  expect_false(is.null(coef(Y)))
  expect_false(any(is.na(fitted(Y))))
  expect_false(any(is.na(resid(Y))))
  expect_false(any(is.na(coef(Y))))
  expect_output(print(Y))
  expect_output(print(summary(Y)))
}

test.wilmoth(W0)
test.wilmoth(W1)

# Build life tables with various choices of 2 input parameters ---
L1 <- wilmothLT(W0, q0_5 = 0.05, k = 0.1)
L2 <- wilmothLT(W0, q0_5 = 0.05, e0 = 65)
L3 <- wilmothLT(W0, q0_5 = 0.05, q15_45 = 0.2)
L4 <- wilmothLT(W0, q0_5 = 0.05, q15_35 = 0.125)
L5 <- wilmothLT(W0, q0_1 = 0.01, k = 0.1)
L6 <- wilmothLT(W0, q0_1 = 0.01, e0 = 65)
L7 <- wilmothLT(W0, q0_1 = 0.05, q15_45 = 0.2)
L8 <- wilmothLT(W0, q0_1 = 0.05, q15_35 = 0.125)
L9 <- wilmothLT(W0, k = 0.01, e0 = 65)
L10 <- wilmothLT(W0, k = 0.01, q15_45 = 0.2)
L11 <- wilmothLT(W0, k = 0.01, q15_35 = 0.125)
L12 <- wilmothLT(W0, q15_45 = 0.125, e0 = 65)
L13 <- wilmothLT(W0, q15_35 = 0.125, e0 = 65)


test.wilmothLT <- function(P){
  expect_true(class(P) == "wilmothLT")
  expect_false(is.null(P))
  expect_false(any(is.na(P$lt)))
  expect_true(all(P$lt[,-1] >= 0))
  expect_false(any(is.na(P$values)))
  
}

for (i in 1:13) test.wilmothLT(get(paste0("L",i)))

# Test Q5 ---
Q5 = 1 - (1 - L1$lt$qx[1])*(1 - L1$lt$qx[2])
expect_equal(Q5, 0.05)


# ----------------------------------------------------------------------------
# Test messages


expect_warning(
  wilmothLT(W0, q15_35 = 0.125, e0 = 65, maxit = 1)
)

expect_error(
  wilmothLT(W0, q15_35 = 0.125, e0 = 65, k = 0.1)
)

expect_error(
  wilmothLT(W0, q15_35 = 0.125, q15_45 = 0.15)
)

expect_error(
  wilmothLT(W0, q0_1 = 0.01, q0_5 = 0.05)
)

expect_error(
  wilmoth(as.character(x), LT = H0)
)

expect_error(
  wilmoth(x[-1], mx = H1)
)

expect_error(
  wilmoth(x[-1], mx = H0)
)

expect_error(
  wilmoth(x, LT = H0, control = list(tol.biweight = -1))
)









