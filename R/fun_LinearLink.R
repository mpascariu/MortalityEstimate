# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# Last update: Fri Apr 30 14:19:06 2021
# --------------------------------------------------- #

#' Fit the Linear-Link Model
#' 
#' Calibrate the Linear-Link mortality model introduced in
#' \insertCite{pascariu2020;textual}{MortalityEstimate}.
#' 
#' @param x Numerical vector containing ages corresponding to the input 
#' data (mx).
#' @param mx Death rates matrix with age as row and time as column.
#' @param y Vector of years corresponding to the mx matrix.
#' @param country Optional. The name of the country that the data 
#' corresponds to. The name is adopted in the output tables.
#' @param theta Age to be fitted.
#' @param use.smooth Logical variable indicating whether the spline smoothing 
#' is applied or not to the estimated coefficients (bx and vx). The smoothing 
#' can be applied in order to avoid jumps in the mortality rates from one 
#' age to another. This using splines. One degree of freedom is allocated 
#' for every 5 year of age.
#' @param method Optimizing method. Least squared approach \code{LSE} or 
#' Poisson likelihood estimation \code{MLE} based on the approach described in 
#' \insertCite{brouhns2002;textual}{MortalityEstimate}. Default: \code{LSE}.
#' @return A \code{LinearLink} object containing:
#'  \item{input}{List with input objects provided in the function}
#'  \item{coefficients}{Estimated coefficient}
#'  \item{fitted.values}{Fitted values}
#'  \item{residuals}{Estimated deviance residuals}
#'  \item{fitted.life.tables}{ Life tables constructed using the fitted values}
#'  \item{df_splines}{Degrees of freedom used in spline smoothing procedure}
#'  \item{model_info}{Description of the model}
#'  \item{process_date}{Data and time stamp}
#' @export
#' @references \insertAllCited{}
#' 
#' Pascariu MD (2018). PhD Thesis: Modelling and Forecasting Mortality.
#' University of Southern Denmark, 53-70. URL:
#' \url{https://github.com/mpascariu/PhD-Thesis/blob/master/Thesis.pdf}.
#' @examples 
#' # Select the 1965 - 1990 time interval and fit the Linear-Link model
#' ages  <- 0:100 # available ages in our datasets
#' years <- 1965:1990 # available years
#' sex   <- 'female'
#' SWEmx <- HMD4mx$SWE[paste(ages), paste(years)]
#' 
#' # Fit the Linear-Link using the least square approach (LSE).
#' M <- LinearLink(x  = ages,
#'                 mx = SWEmx,
#'                 y  = years,
#'                 country = 'SWEDEN',
#'                 theta   = 0,
#'                 method  = 'LSE')
#' M
#' summary(M)
#' coef(M)
#' ls(M)
#' 
#' 
#' # Derive a mortality curve (life table) from a value of
#' # life expectancy at birth in 2014, say 84.05
#' e0 <- 84.05
#' LT1 <- LinearLinkLT(M, ex = e0)
#' LT2 <- LinearLinkLT(M, ex = e0, use.vx.rotation = TRUE) 
LinearLink <- function(x, 
                       mx, 
                       y,
                       country = '...', 
                       theta = 0,
                       use.smooth = TRUE, 
                       method = 'LSE'){
  #-------------------------------------------------
  input <- c(as.list(environment()))
  check_LinearLink_input(input)
  
  cat('\n   Fitting Linear-Link model\n')
  pb <- startpb(0, length(y)) # Start the clock!
  on.exit(closepb(pb))        # Stop clock on exit.
  #-------------------------------------------------
  # Data preparation
  mx_input   <- as.matrix(mx)
  dimnames(mx_input) <- list(x, y)
  model_info <- "Linear-Link: log m[x] = b[x] log e[x] + k v[x]"
  # Compute multiple life tables (in order to get ex)
  LT <- data.frame()
  for (i in 1:ncol(mx)) {
    LT_i <- LifeTable(x, mx = mx_input[, i])$lt
    LT_i <- cbind(country = country, year = y[i], LT_i,
                  ex0 = LT_i[LT_i$x == theta, 'ex'], Ex = 1)
    LT_i <- LT_i[complete.cases(LT_i), ]
    LT   <- rbind(LT, LT_i)
  }
  #-------------------------------------------------
  # Step 1   - Takes place before entering this function.
  # Step 2-3 - Estimate bx and vx
  log_ex_theta <- log(LT[LT$x == theta, 'ex'])
  log_mx       <- t(log(mx_input))
  if (method == 'LSE') fit_link <- fitw_LSE(log_ex_theta, log_mx)
  if (method == 'MLE') fit_link <- fitw_MLE(log_ex_theta, log_mx)
  k_ <- fit_link$k_
  #-------------------------------------------------
  # Step 4 - Smooth Coefficients
  coeffs_raw    <- data.frame(bx = fit_link$bx, vx = fit_link$vx, row.names = x)
  coeffs        <- coeffs_raw
  coeffs_smooth <- coeffs_raw*0
  degrees       <- ifelse(use.smooth, round(length(x)/5), x)
  df_spline     <- ifelse(use.smooth, degrees, 'Smoothing not used')
  
  if (use.smooth) {
    for (j in 1:ncol(coeffs_raw)) {
      coeffs_smooth[, j] <- smooth.spline(coeffs_raw[, j], df = degrees)$y
    }
    coeffs_smooth[1, 1] <- coeffs_raw[1, 1] # leave infant mortality unsmoothed
    coeffs <- coeffs_smooth
  }
  #--------------------------------------------------
  # Step 5-6 - Compute fitted values of mx using precise k
  tab_ex   <- LT[LT$x == theta, c('year', 'ex')]
  LT_optim <- NULL
  for (i in 1:length(y)) {
    optim_obj  <- compute_lt_optim(x, coeffs, tab_ex[i, 2])
    LT_optim_i <- cbind(country = country, year = tab_ex[i, 1], optim_obj$lt)
    LT_optim   <- rbind(LT_optim, LT_optim_i)
    k_[i]      <- optim_obj$k
    setpb(pb, i)
  }
  fitted_mx <- reshape(
    data = LT_optim[, c("country", "year", "x", "mx")], 
    direction = 'wide', 
    idvar = c('country', 'x'), 
    timevar = 'year')[, -(1:2)]
  dimnames(fitted_mx) <- list(x, y)
  residuals    <- mx_input - fitted_mx
  coefficients <- list(bx = coeffs$bx, vx = coeffs$vx, k = k_)
  #-----------------------------------
  # Output
  out <- list(input = input,
              call = match.call(),
              coefficients = coefficients, 
              fitted = fitted_mx,
              residuals = residuals, 
              fitted.life.tables = LT_optim,
              df_spline = df_spline, 
              model_info = model_info, 
              process_date = date())
  out <- structure(class = 'LinearLink', out)
  return(out)
}

#' Check consistency in input arguments
#' @param input a list containing the input objects of the LinearLink function
#' @keywords internal
check_LinearLink_input <- function(input){
  with(input, {
    
    if (nrow(mx) != length(x)) {
      stop("\nnrow(mx) must be equal to length(x)", call. = FALSE)
    }
    
    if (ncol(mx) != length(y)) {
      stop("\nncol(mx) must be equal to length(y)", call. = FALSE)
    }
    
    if (theta > 50 & method == 'LSE') {
      message(
        "For theta > 50 the MLE method has been observed to be more reliable."
      )
    }
    if (!(method %in% c('LSE', 'MLE'))) {
      stop(paste("Method", method, "not available. Try 'LSE' or 'MLE' "), 
           call. = FALSE)
    }
  })
}



#' Function to optimize a life table
#' 
#' @param ages ages
#' @param coefs coefficients
#' @param ex0 life expectancy
#' @keywords internal
compute_lt_optim <- function(x, coefs, ex0){
  penalty <- function(k){
    mx.hat  <- exp(coefs[, 1] * log(ex0) + coefs[, 2] * k)
    ex0.hat <- LifeTable(x = x, mx = mx.hat)$lt$ex[1]
    abs(ex0.hat - ex0)
  }
  
  k.hat  <- optim(0, penalty, method = "Brent", upper = 150, lower = -250)$par
  mx.hat <- exp(coefs[, 1]*log(ex0) + coefs[, 2]*k.hat)
  lt.hat <- LifeTable(x = x, mx = mx.hat)$lt
  out    <- list(k = k.hat, lt = lt.hat)
  return(out)
}


# Two fitting procedures for Linear-Link model
# 1. LSE + SVD
# 2. Poisson MLE

# 1. ------------------------------------------------------
#' Estimate bx using least square method and vx with SVD
#' 
#' @param log_ex_theta Life expectancy at age theta (log values)
#' @param log_mx death rates (log values)
#' @inheritParams wilmoth_control
#' @keywords internal
#'
fitw_LSE <- function(log_ex_theta, log_mx, nu = 1, nv = 1){
  # Step 2 - Fit bx
  fit <- LSEfit(x = log_ex_theta, y = log_mx)
  bx  <- as.numeric(fit$coef)
  # Step 3 - Compute residuals and fit SVD portion of model
  fitted_log_mx <- log_ex_theta %*% t(bx)
  resid_log_mx  <- fitted_log_mx - log_mx
  dimnames(fitted_log_mx) <- dimnames(resid_log_mx) <- dimnames(log_mx)
  resid_log_mx[resid_log_mx == Inf] <- unique(sort(x = resid_log_mx, 
                                                   decreasing = TRUE))[2]
  vx  <- svd(resid_log_mx, nu, nv)$v
  if (min(vx) < 0) { 
    # Shift upwards the vx curve if negative values found
    vx <- vx + abs(min(vx)) 
  }
  vx  <- vx / sum(vx) # scale to 1
  k_  <- rep(NA, length(log_ex_theta))
  
  #Exit
  out <- list(bx = bx, 
              vx = vx, 
              k_ = k_)
  return(out)
}


#' Find the Least Squares Fit
#' 
#' @inheritParams bifit
#' @inheritParams wilmoth_control
#' @keywords internal
LSEfit <- function(x, 
                   y, 
                   intercept = FALSE, 
                   tol.lsfit = 1e-07) {
  
  if (is.vector(y)) {
    z   <- lsfit(x = x, 
                 y = y, 
                 wt = NULL, 
                 intercept = intercept, 
                 tolerance = tol.lsfit)
    coef  <- z$coef
    resid <- z$resid
  }
  
  if (is.matrix(y)) {
    resid <- coef <- NULL
    for (j in 1:ncol(y)) {
      Y     <- y[, j] 
      Y[Y == -Inf] <- -10
      z     <- LSEfit(x, Y, intercept, tol.lsfit)
      resid <- cbind(resid, z$resid)
      coef  <- cbind(coef, z$coef) 
    }
  }
  
  out <- list(coef = coef, residuals = resid)
  return(out) 
}


# 2. ------------------------------------------------------
#' # Estimate bx, vx and k with Poisson Likelihood Method 
#' 
#' The implemented method of estimated the bx, vx and k coefficients of the 
#' LinearLink model is based on the approach described in 
#' \insertCite{brouhns2002;textual}{MortalityEstimate}
#' for fitting the Lee-Carter model. 
#' Code written by Jose Manuel Aburto with minor changes by Marius Pascariu
#' @references \insertAllCited{}
#' @keywords internal
fitw_MLE <- function(log_ex_theta, log_mx, ...){
  # Normally, deaths and exposes is needed in order to fit the model using 
  # the Poisson distribution. However, if a vector of mx is available we can 
  # estimate Dx (deaths) and Ex (exposures) in such a way that the parameters 
  # are reasonable computed.
  Dx  <- t(exp(log_mx)) * 1e6 # Dx estimation
  Ex  <- Dx*0 + 1e6           # Ex estimation
  fit <- PoissonMLE(log_ex_theta, Dx, Ex, ...)
  bx  <- fit$bx
  vx  <- matrix(fit$vx, ncol = 1)
  if (min(vx) < 0) vx <- vx + abs(min(vx)) # Shift up vx curve if negative
  k_  <- as.numeric(fit$k)
  out <- list(bx = bx, vx = vx, k_ = k_)
  return(out)
}


#' @keywords internal
#'
PoissonMLE <- function(log_ex_theta, 
                       Dx, 
                       Ex, 
                       iter = 500, 
                       tol = 1e-04){
  
  # dimensions
  n <- ncol(Dx)
  # Initialise
  mat_1    <- matrix(1, nrow = n, ncol = 1)    
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
    temp   <- update_alpha(alpha, vx, k, Dx, Ex, Dx_fit, mat_1)
    Dx_fit <- temp$Dx_fit
    alpha  <- temp$alpha
    #
    temp   <- update_vx(alpha, vx, k, Dx, Ex, Dx_fit, mat_1)
    Dx_fit <- temp$Dx_fit
    vx     <- temp$vx
    #
    temp   <- update_k(alpha, vx, k, Dx, Ex, Dx_fit, mat_1)
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
  
  vxk <- vx %*% t(k)
  log.MU.hat <- alpha %*% t(mat_1) + vxk
  bx_hat <- rowMeans((log.MU.hat - vxk)/log_ex_theta)
  
  # output
  out <- list(bx = as.numeric(bx_hat), 
              vx = as.numeric(vx), 
              k = as.numeric(k))
  return(out)
}

#' Update alpha
#' @keywords internal
update_alpha <- function(alpha, vx, k, Dx, Ex, Dx_fit, mat_1){
  difD   <- Dx - Dx_fit
  alpha  <- alpha + difD %*% mat_1 / (Dx_fit %*% mat_1)
  Eta    <- alpha %*% t(mat_1) + vx %*% t(k)
  Dx_fit <- Ex * exp(Eta)
  list(alpha = alpha, Dx_fit = Dx_fit)
}

#' Update vx
#' @keywords internal
update_vx <- function(alpha, vx, k, Dx, Ex, Dx_fit, mat_1){
  difD   <- Dx - Dx_fit  # exp(log_mx) - Dx_fit
  vx     <- vx + difD %*% k / (Dx_fit %*% (k ^ 2))
  Eta    <- alpha %*% t(mat_1) + vx %*% t(k)
  Dx_fit <- Ex * exp(Eta)
  list(vx = vx, Dx_fit = Dx_fit)
}

#' Update k
#' @keywords internal
update_k <- function(alpha, vx, k, Dx, Ex, Dx_fit, mat_1){
  difD   <- Dx - Dx_fit
  k      <- k + t(difD) %*% vx / (t(Dx_fit) %*% (vx ^ 2))
  k      <- k - mean(k)
  k      <- k / sqrt(sum(k ^ 2))
  k      <- matrix(k, ncol = 1)
  Eta    <- alpha %*% t(mat_1) + vx %*% t(k)
  Dx_fit <- Ex * exp(Eta)
  list(k = k, Dx_fit = Dx_fit)
}



#' Estimate LinearLink Life Table
#' 
#' Construct a life table based on the Linear-Link estimates and a given value 
#' of life expectancy at age theta 
#' (\insertCite{pascariu2020;textual}{MortalityEstimate}).
#' @param object An object of class \code{LinearLink}.
#' @param ex Value of life expectancy for which we want to estimate 
#' the mortality curve. Type: numerical scalar.
#' @param use.vx.rotation Logical argument. If \code{TRUE} the adjustment method
#' inspired from \insertCite{li2013;textual}{MortalityEstimate} article is 
#' applied to the vx coefficients before estimated the life table. 
#' If \code{FALSE} the fitted vx coefficients are used in the estimations of 
#' the life table.
#' @param ... Additional arguments affecting the predictions produced.
#' @return Predicted values of the Linear-Link model.
#' @references \insertAllCited{}
#' @seealso \code{\link{LinearLink}}
#' @examples 
#' # See examples in the ?LinearLink help page
#' @export
LinearLinkLT <- function(object, 
                         ex, 
                         use.vx.rotation = FALSE, 
                         ...) {
  
  # Choose vx coefficients
  x  <- object$input$x
  bx <- coef(object)$bx
  
  if (use.vx.rotation == TRUE) {
    
    if (object$input$theta > 0) {
      stop("Currently vx rotation is implemented only for theta = 0.",
           " Set use.vx.totation = FALSE.", call. = FALSE)
    }
    vx = rotate_vx(object, ex_target = ex, ...) 
    
  } else { 
    vx = coef(object)$vx 
  }
  # Data.frame with all coefficients used in prediction
  coefs <- data.frame(bx = bx, vx = vx, row.names = x)
  # Find the right life table
  out <- compute_lt_optim(x, coefs = coefs, ex0 = ex)
  out$bx <- bx
  out$vx <- vx
  
  return(out)
}


#' Compute Rotated vx Coefficients
#' 
#' This functions computes the rotation of the vx coefficients using the method
#' presented in \insertCite{li2013;textual}{MortalityEstimate} paper. 
#' @param object An object of class 'LinearLink'
#' @param ex_target A value of life expectancy for which we want to derive 
#' the rotated vx
#' @param e0_u Ultimate value of life expectancy. At this point the rotation 
#' process reaches its maximum efficiency. Here. e0_u = 80 is taken as the 
#' default.  
#' @param e0_threshold Level of life expectancy where the rotation should begin.
#' If rotated_vx is computed for ex_target <= e0_threshold then no difference
#' will be observed.
#' @param e0_conv Set limit for convergence. For life expectancy a birth 
#' the default value is 130.  
#' @param p_ The power to the smooth-weight function, p_, takes values 
#' between 0 and 1, which makes the rotation faster at starting times and 
#' slower at ending times. Here, p = .5 is taken as the default
#' @references \insertAllCited{}
#' @keywords internal
#' 
rotate_vx <- function(object, 
                      ex_target, 
                      e0_u = 102, 
                      e0_threshold = 80,
                      e0_conv = 130,
                      p_ = 0.5){
  
  vx = coef(object)$vx
  x  = object$input$x
  x1 = max(min(x), 15):65 # young ages
  x2 = 66:max(x) # old ages
  
  vx_u = vx * 0
  vx_young_ages = mean(vx[x1 + 1]) # select vx corresponding to young ages. 
  # Only the values between age 15 and 65 are being used.
  
  # Derive a logistic shape. The values have to be between 0 and 1, 
  # they will be scaled.
  n_vx        <- length(vx[x2]) # count age groups in x2
  n_vx_ext    <- n_vx + (e0_conv - max(x)) # set limit for convergence at 130
  x_num       <- seq(-6, 6, length.out = n_vx_ext) 
  logit_shape <- 1 - exp(x_num) / (1 + exp(x_num)) 
  vx_old_ages <- logit_shape * vx_young_ages # scale values
  # This is our vx ultimate (vx_u)
  vx_u[min(x):max(x1 + 1)] <- vx_young_ages
  vx_u[x2 + 1] <- vx_old_ages[1:n_vx]
  vx_u <- vx_u/sum(vx_u)  ## rescale
  # Compute weights
  w_t  <- (ex_target - e0_threshold)/(e0_u - e0_threshold)
  ws_t <- (0.5 * (1 + sin(pi/2 * (2*w_t - 1))) ) ^ p_ 
  # Compute rotate_vx
  if ( ex_target < e0_threshold ) rot_vx <- vx 
  if ( e0_threshold <= ex_target & ex_target < e0_u ) {
    rot_vx <- (1 - ws_t) * vx + ws_t * vx_u
  }
  if (ex_target >= e0_u) rot_vx <- vx_u 
  return(rot_vx)
}


# S Functions
# Classes and methods

#' @keywords internal
#' @export
summary.LinearLink <- function(object, ...) {
  mi    <- object$model_info
  cl    <- object$call
  dev   <- round(summary(as.vector(as.matrix(object$residuals))), 5)
  coefs <- data.frame(bx = coef(object)$bx, 
                      vx = coef(object)$vx,
                      row.names = object$input$x)
  Hbxvx <- head_tail(coefs, digits = 5, hlength = 6, tlength = 6)
  Hk    <- head_tail(data.frame(. = '.', k = coef(object)$k),
                     digits = 5, 
                     hlength = 6, 
                     tlength = 6)
  H   <- data.frame(Hbxvx, Hk)
  dfs <- object$df_spline
  
  out <- list(model_info = mi, 
              call = cl, 
              dev = dev, 
              coef = H, 
              df_spline = dfs)
  out <- structure(class = "summary.LinearLink", out)
  return(out)
}

#' @keywords internal
#' @export
print.summary.LinearLink <- function(x, ...) {
  cat('\nModel:\n'); cat(x$model_info)
  cat("\n\nCall:\n"); print(x$call)
  cat('\nDeviance Residuals:\n'); print(x$dev)
  cat("\nCoefficients:\n"); print(x$coef)
  cat("\nSpline Smoothing (degrees of freedom): "); cat(x$df_spline)
}

#' @keywords internal
#' @export
print.LinearLink <- function(x, ...){
  cat(x$model_info, "\n")
  with(x$input,
       {
         cat('\nFitted for life expectancy at age:', theta)
         cat('\nTime interval:', min(y), '-', max(y))
         cat('\nAge-range:', min(x), '-', max(x))
         cat('\nCountry:', country, '\n')
         met <- ifelse(method == 'LSE', 'Least Squares (LSE)', 
                       'Poisson Maximum Likelihood (MLE)')
         cat('\nFitting Procedure:', met)
         cat('\nSmoothing:', use.smooth)
       })
}
