# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Dec  4 22:16:55 2018
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

W1$k

test.wilmoth <- function(Y) {
  expect_true(class(Y) == "wilmoth")
  expect_false(is.null(Y))
  expect_false(is.null(coef(Y)))
  expect_false(any(is.na(fitted(Y))))
  expect_false(any(is.na(resid(Y))))
  expect_false(any(is.na(coef(Y))))
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
