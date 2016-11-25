###############  Fit linear model using least square method  ###################
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

#' @keywords internal
#'
FUN.mxhat <- function(ages, coefs, ex0 = 0, k = 0) {
     mx <- exp(coefs[, 1]*log(ex0) + coefs[, 2]*k)
     return(mx) 
}

#' @keywords internal
#'
FUN.lt_k0 <- function(ages, coefs, ex0 = 0, k = 0) {
     fx <- FUN.mxhat(ages, coefs, ex0, k)
     LT <- lifetable(ages, mx = fx)
     lt <- LT$lt
     lt.exact <- LT$lt.exact
     out <- list(lt = lt, lt.exact = lt.exact)
     return(out)
}


#' Function to optimize a life table
#' 
#' @param ages ages
#' @param coefs coefficients
#' @param ex0 life expectancy
#' @return Results
#' @keywords internal
#' 
FUN.lt_optim <- function(ages, coefs, ex0){
     ptm <- proc.time() # Start the clock!
     penalty <- function(k_init){
          ex_k <- FUN.lt_k0(ages, coefs, ex0, k = k_init)$lt.exact$ex[1]
          out  <- abs(ex_k - ex0)
          return(out)
     }
     k.optim <- optim(0, penalty, method = "Brent", 
                      upper = 30, lower = -30)$par
     # LT.init <- FUN.lt_k0(ages, coefs, ex0, k = 0)
     LT      <- FUN.lt_k0(ages, coefs, ex0, k = k.optim)
     proc_speed <- round((proc.time() - ptm)[3], 2) # Stop the clock
     out <- list(k = k.optim, lt = LT$lt, lt.exact = LT$lt.exact, 
                 process_date = date(), process_speed = proc_speed)
     return(out)
}

#' Fit Linear-Link model
#' 
#' @param mx Death rates matrix with age as row and time as column
#' @param mx_ages vector of ages corresponding to the mx matrix
#' @param mx_years vector of years corresponding to the mx matrix
#' @param mx_country The name of the country that the data corresponds to.
#' The name will be adopted in the output tables. It is optional.
#' @param x_fit Age to be fitted
#' @param use.smooth Logical variable indicating wheter the spline smoothing is 
#' applyed or not to the estimated coefficients (bx and vx).
#' @param method Optimizing method. Least squared approch \code{LSE} or Poisson 
#' likelihood estimation \code{MLE}.
#' @return A \code{LinearLink} object
#' @export
#' @examples 
#' library(LinearLink)
#'  
#' # Select the 1965 - 1990 time interval and fit the Linear-Link model
#' ages    <- 0:100 # available ages in our datasets
#' years   <- 1965:1990 # available years
#' 
#' SWEmx <- HMD.test.data$SWE[paste(ages), paste(years)]
#' 
#' fit_SWE <- LinearLink(mx = SWEmx, mx_ages = ages, mx_years = years, 
#'                       mx_country = 'SWEDEN', x_fit = 0, method = 'LSE')
#' summary(fit_SWE)
#' 
#' # Derive a mortality curve (life table) from a value of
#' # life expectancy at birth in 2014, say 84.05
#' new_e0   <- 84.05
#' pred_SWE <- predict(fit_SWE, new_e0)
#' pred_SWE$lt
#' 
LinearLink <- function(mx, mx_ages, mx_years, mx_country = '...', 
                           x_fit = 0, use.smooth = TRUE, method = 'LSE'){
  input <- c(as.list(environment()))
  pb <- startpb(0, length(mx_years)) # Start the clock!
  on.exit(closepb(pb)) # Stop clock on exit.
  
  # Data preparation ---------------------------------
  model_name <- "Linear-Link (2016): ln[m(x)] = b(x)ln[e(x)] + kv(x)"
  mx_input <- as.matrix(mx)
  dimnames(mx_input) <- list(mx_ages, mx_years)
  
  # Compute life expectancy --------------------------
  LT <- data.frame()
  for (i in 1:ncol(mx)) {
    year_i <- mx_years[i]
    LT_i   <- lifetable(x = mx_ages, mx = mx_input[, i])$lt
    LT_i   <- cbind(country = mx_country, year = year_i, LT_i, 
                    ex0 = LT_i[LT_i$age == x_fit, 'ex'], Ex = 1)
    LT_i   <- LT_i[complete.cases(LT_i), ]
    LT <- rbind(LT, LT_i)
  }
  #-------------------------------------------------
  # Fit linear portion of model
  log_ex_theta <- log(LT[LT$age == x_fit, 'ex'])
  log_mx_theta <- t(log(mx_input))
  if (method == 'LSE') {
    fit_link.LSE <- FUN.bifit(y = log_mx_theta, x = log_ex_theta)
    bx <- as.numeric(fit_link.LSE$coef)
  }
  if (method == 'MLE') {
    fit_link.MLE <- suppressWarnings( 
      gnm(mx ~ -1 + as.factor(age):log(ex0) + offset(log(Ex)), 
          family = poisson(link = "log"), data = LT))
    bx <- as.numeric(coef(fit_link.MLE))
  }
  #-------------------------------------------------
  # Compute residuals and fit SVD portion of model
  fitted_log_mx <- log_ex_theta %*% t(bx)
  dimnames(fitted_log_mx) <- dimnames(log_mx_theta)
  resid_log_mx <- fitted_log_mx - log_mx_theta
  resid_log_mx[resid_log_mx == Inf] <- unique(sort(resid_log_mx, decreasing = T))[2]
  vx  <- svd(resid_log_mx, 1, 1)$v
  
  # Coefficients ----------------------------------
  coefs.raw <- round(data.frame(bx , vx, row.names = mx_ages), 8)
  coeffs <- coefs.raw
  # I have to smooth the vx's so that we can avoid jumps 
  # in the mortality rates from one age to another. We can do this
  # by using splines. We can allocate 1 degree of freedom for every 5 ages.
  degrees   <- ifelse(use.smooth, round(length(mx_ages)/5), mx_ages)
  df_spline <- ifelse(use.smooth, degrees, 'Smoothing not used')
  
  if (use.smooth) {
    coefs.smooth <- coefs.raw*0
    for (j in 1:ncol(coefs.raw)) {
      coefs.smooth[, j] <- smooth.spline(coefs.raw[, j], df = degrees)$y
    }
    coefs.smooth[1, 1] <- coefs.raw[1, 1] # leave infant mortality unsmoothed.
    coeffs <- coefs.smooth
  }
  #--------------------------------------------------
  # Compute fitted values of mx ---- 
  table_ex <- LT[LT$age == x_fit, c('year', 'ex')]
  k_values = LT_optim <- NULL
  for (i in 1:length(mx_years)) {
    year_i <- table_ex[i, 1]
    ex_target_i <- table_ex[i, 2]
    Optim_out <- FUN.lt_optim(ages = mx_ages, coefs = coeffs, 
                              ex0 = ex_target_i)
    LT_optim_i <- cbind(country = mx_country, 
                        year = year_i, Optim_out$lt)
    colnames(LT_optim_i)[3] <- 'age'
    LT_optim <- rbind(LT_optim, LT_optim_i)
    k_optim  <- data.frame(country = mx_country, year = year_i, 
                           ex = LT_optim_i[LT_optim_i$age == x_fit, ]$ex, 
                           k = round(Optim_out$k, 6))
    k_values <- rbind(k_values, k_optim)
    setpb(pb, i)
  }
  
  table_mx <- LT_optim[, 1:4]
  fitted.values <- spread(table_mx, key = year, value = mx)[, -(1:2)]
  rownames(fitted.values) <- mx_ages
  residuals <- mx_input - fitted.values
  coefficients <- list(bx = coeffs$bx, vx = coeffs$vx, k = k_values$k)
  #-----------------------------------
  # Output
  out = structure(class = 'LinearLink',
                 list(input = input, df_spline = df_spline, 
                      coefficients = coefficients, fitted.values = fitted.values, 
                      residuals = residuals, fitted.life.tables = LT_optim, 
                      model_name = model_name, process_date = date()))
  out$call <- match.call()
  return(out)
}









