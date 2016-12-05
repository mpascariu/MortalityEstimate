# This file 2 fitting procedures for Linear-Link model
# 1. LSE + SVD
# 2. Poisson MLE


# 1. ------------------------------------------------------
#' Estimate bx using least square method and vx with SVD
#' @keywords internal
#'
fitw_LSE <- function(log_ex_theta, log_mx, ...){
  # Step 2 - Fit bx
  fit = FUN.bifit(x = log_ex_theta, y = log_mx)
  bx  = as.numeric(fit$coef)
  # Step 3 - Compute residuals and fit SVD portion of model
  fitted_log_mx <- log_ex_theta %*% t(bx)
  dimnames(fitted_log_mx) <- dimnames(log_mx)
  resid_log_mx <- fitted_log_mx - log_mx
  resid_log_mx[resid_log_mx == Inf] <- unique(sort(resid_log_mx, 
                                                   decreasing = T))[2]
  vx  <- svd(resid_log_mx, 1, 1)$v
  vx  <- vx / sum(vx) # scale to 1
  k_  <- rep(NA, length(log_ex_theta))
  out <- list(bx = bx, vx = vx, k_ = k_)
  return(out)
}

#'@keywords internal
#'
FUN.bifit <- function(x, y) {
  
  if (is.vector(y)) {
    z <- lsfit(x, y, intercept = FALSE)
    z.resid  <- z$resid
    coef.new <- z$coef
    return(list(coef = coef.new, residuals = z.resid)) }
  
  if (is.matrix(y)) {
    resid <- coef <- NULL
    for (j in 1:ncol(y)) {
      y.j   <- y[, j] 
      y.j[y.j == -Inf] = -10
      z     <- FUN.bifit(x, y.j)
      resid <- cbind(resid, z$resid)
      coef  <- cbind(coef, z$coef) 
    }
    out <- list(coef = coef, residuals = resid)
    return(out) 
  }
}


# 2. ------------------------------------------------------
#' # Estimate bx, vx and k with Poisson Likelihood Method 
#' 
#' The implemeted method of estimated the bx, vx and k coefficients of the 
#' LinearLink model is based on the approach described in Brouhns et al. 2002
#' for fittin the Lee-Carter model. 
#' Code writen by Jose Manuel Aburto with minor changes by Marius Pascariu
#' @source Brouhns et al. 2002
#' @keywords internal
#'
fitw_MLE <- function(log_ex_theta, log_mx, ...){
  # Normally, deaths and exposes is needed in order to fit the model using 
  # the Poisson distribution. However if a vector of mx is available we can 
  # estimate Dx (deaths) and Ex (exposures) in such a way that the parameters 
  # are resonable computed.
  Dx = t(exp(log_mx)) * 1e6 # Dx estimation
  Ex = Dx*0 + 1e6 # Ex estimation
  fit <- PoissonMLE(log_ex_theta, Dx, Ex, ...)
  out <- list(bx = fit$bx, vx = matrix(fit$vx, ncol = 1), 
              k_ = as.numeric(fit$k))
  return(out)
}


#' @keywords internal
#'
PoissonMLE <- function(log_ex_theta, Dx, Ex, iter = 500, tol = 1e-04){
  # dimensions
  n <- ncol(Dx)
  # Initialise
  mat_1    <- matrix(1, nrow = ncol(Dx), ncol = 1)    
  Fit.init <- log((Dx + 1)/(Ex + 2))
  Dx_fit   <- Ex * exp(Fit.init)  # Ex * exp(log_mx)
  alpha    <- Fit.init %*% mat_1 / n
  
  vx     <- matrix(1 * alpha, ncol = 1)
  sum_vx <- sum(vx) 
  vx     <- vx / sum_vx
  
  k <- matrix(seq(n, 1, by = -1), nrow = n, ncol = 1)
  k <- k - mean(k)
  k <- k / sqrt(sum(k * k))
  k <- k * sum_vx
  
  # Iteration
  for (i in 1:iter) {
    alpha_old = alpha; vx_old = vx; k_old = k
    #
    temp   <- Update.alpha(alpha, vx, k, Dx, Ex, Dx_fit)
    Dx_fit <- temp$Dx_fit
    alpha  <- temp$alpha
    #
    temp   <- Update.vx(alpha, vx, k, Dx, Ex, Dx_fit)
    Dx_fit <- temp$Dx_fit
    vx     <- temp$vx
    #
    temp   <- Update.k(alpha, vx, k, Dx, Ex, Dx_fit)
    Dx_fit <- temp$Dx_fit
    k      <- temp$k
    crit   <- max(max(abs(alpha - alpha_old)), 
                  max(abs(vx - vx_old)), 
                  max(abs(k - k_old)))
    if (crit <= tol) break
  }
  
  #constraints
  k  <- k * sum(vx)
  vx <- vx / sum(vx) # scale to 1

  log.MU.hat <- alpha %*% t(mat_1) + vx %*% t(k)
  bx_hat <- rowMeans((log.MU.hat - vx %*% t(k))/log_ex_theta)
  
  # output
  out <- list(bx = as.numeric(bx_hat), 
              vx = as.numeric(vx), k = as.numeric(k))
  return(out)
}

#' Update alpha
#' @keywords internal
#'
Update.alpha <- function(alpha, vx, k, Dx, Ex, Dx_fit){
  mat_1 <- matrix(1, nrow = ncol(Dx), ncol = 1) 
  difD  <- Dx - Dx_fit
  alpha <- alpha + difD %*% mat_1 / (Dx_fit %*% mat_1)
  Eta   <- alpha %*% t(mat_1) + vx %*% t(k)
  Dx_fit <- Ex * exp(Eta)
  list(alpha = alpha, Dx_fit = Dx_fit)
}

#' Update vx
#' @keywords internal
#'
Update.vx <- function(alpha, vx, k, Dx, Ex, Dx_fit){
  mat_1 <- matrix(1, nrow = ncol(Dx), ncol = 1) 
  difD  <- Dx - Dx_fit  # exp(log_mx) - Dx_fit
  vx    <- vx + difD %*% k / (Dx_fit %*% (k ^ 2))
  Eta   <- alpha %*% t(mat_1) + vx %*% t(k)
  Dx_fit <- Ex * exp(Eta)
  list(vx = vx, Dx_fit = Dx_fit)
}

#' Update k
#' @keywords internal
#'
Update.k <- function(alpha, vx, k, Dx, Ex, Dx_fit){
  mat_1 <- matrix(1, nrow = ncol(Dx), ncol = 1) 
  difD <- Dx - Dx_fit
  k <- k + t(difD) %*% vx / (t(Dx_fit) %*% (vx ^ 2))
  k <- k - mean(k)
  k <- k / sqrt(sum(k ^ 2))
  k <- matrix(k, ncol = 1)
  Eta <- alpha %*% t(mat_1) + vx %*% t(k)
  Dx_fit <- Ex * exp(Eta)
  list(k = k, Dx_fit = Dx_fit)
}
