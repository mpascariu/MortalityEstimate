
#' Fit model using least square method
#' @keywords internal
#'
fitw_LSE <- function(log_ex_theta, log_mx){
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
  vx <- vx / sum(vx) # scale to 1
  
  out <- list(bx = bx, vx = vx)
  return(out)
}

#'@keywords internal
#'
FUN.bifit <- function(x, y) {
  
  if (is.vector(y)) {
    z <- lsfit(x, y, intercept = FALSE)
    z.resid <- z$resid
    coef.new <- z$coef
    return(list(coef = coef.new, residuals = z.resid)) }
  
  if (is.matrix(y)) {
    resid <- coef <- NULL
    for (j in 1:ncol(y)) {
      y.j <- y[, j] 
      y.j[y.j == -Inf] = -10
      z <- FUN.bifit(x, y.j)
      resid <- cbind(resid, z$resid)
      coef <- cbind(coef, z$coef) 
    }
    return(list(coef = coef, residuals = resid)) 
  }
}


# Fit model with Poisson Likelihood Method --------------
#' @keywords internal
#'
fitw_MLE <- function(log_ex_theta, log_mx, LT){
  # Step 2 - Fit bx
    fit <- suppressWarnings( 
      gnm(LT$mx ~ -1 + as.factor(LT$age):log(LT$ex0) + offset(log(LT$Ex)), 
          family = poisson(link = "log")))
    bx <- as.numeric(coef(fit))
    # Step 3 - Compute residuals and fit SVD portion of model
    fitted_log_mx <- log_ex_theta %*% t(bx)
    dimnames(fitted_log_mx) <- dimnames(log_mx)
    resid_log_mx <- fitted_log_mx - log_mx
    resid_log_mx[resid_log_mx == Inf] <- unique(sort(resid_log_mx, 
                                                     decreasing = T))[2]
    vx  <- svd(resid_log_mx, 1, 1)$v
    vx  <- vx / sum(vx) # scale to 1
    
    out <- list(bx = bx, vx = vx)
    return(out)
}

#' @keywords internal
#'
fitw_MLE2 <- function(log_ex_theta, log_mx){
  Dx = t(exp(log_mx)) * 1e6
  Ex = Dx*0 + 1e6
  fit <- PoissonMLE(log_ex_theta, Dx, Ex)
  out <- list(bx = fit$bx, vx = matrix(fit$vx, ncol = 1), k = fit$k)
  return(out)
}



# Based mainly in Brouhns et al 2002
# Some useful functions to fit the model 
# We just reparametrize to fit Pascariu et al's model
# alpha= b_x log(e(theta))
# vx = v_x
# k = k
# Note we need Deaths an Exposures

#' @keywords internal
#'
PoissonMLE <- function(log_ex_theta, Dx, Ex){
  # dimensions
  n <- ncol(Dx)
  # Initialise
  mat_1 <- matrix(1, nrow = ncol(Dx), ncol = 1)    
  Fit.init <- log((Dx + 1)/(Ex + 2))
  Dx_fit <- Ex * exp(Fit.init)  # Ex * exp(log_mx)
  alpha <- Fit.init %*% mat_1 / n
  
  vx <- matrix(1 * alpha, ncol = 1)
  sum_vx <- sum(vx) 
  vx <- vx / sum_vx
  
  k <- matrix(seq(n, 1, by = -1), nrow = n, ncol = 1)
  k <- k - mean(k)
  k <- k / sqrt(sum(k * k))
  k <- k * sum_vx
  # Iteration
  
  for (iter in 1:50) {
    alpha.old <- alpha
    vx.old <- vx
    k.old <- k
    #
    temp <- Update.alpha(alpha, vx, k, Dx, Ex, Dx_fit)
    Dx_fit <- temp$Dx_fit
    alpha <- temp$alpha
    #
    temp <- Update.vx(alpha, vx, k, Dx, Ex, Dx_fit)
    Dx_fit <- temp$Dx_fit
    vx <- temp$vx
    #
    temp <- Update.k(alpha, vx, k, Dx, Ex, Dx_fit)
    Dx_fit <- temp$Dx_fit
    k <- temp$k
    crit <- max(max(abs(alpha - alpha.old)),
                max(abs(vx - vx.old)),
                max(abs(k - k.old)))
    if (crit < 1e-04) break
  }
  
  #constraints
  k  <- k * sum(vx)
  vx <- vx / sum(vx) # scale to 1

  log.MU.hat <- alpha %*% t(mat_1) + vx %*% t(k)
  bx_hat <- rowMeans((log.MU.hat - vx %*% t(k))/log_ex_theta)
  
  # output
  out <- list(bx = as.numeric(bx_hat), vx = as.numeric(vx), k = as.numeric(k))
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
