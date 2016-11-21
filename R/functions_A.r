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


#' Function to create a life table with input variables: (age, Dx, Ex) 
#' or (age, mx) or (age, qx)
#' 
#' @param x Ages
#' @param Dx Death counts
#' @param Ex Exposure
#' @param mx death rates
#' @param qx 1 year death propabilities
#' @param lx0 Radix
#' @param ax0 0.1
#' @return Results
#' @export
#' @examples 
#' library(LinearLink)
#' 
#' mx_vector <- HMD_mx$USA[,1]
#' life.table(x = 0:100, mx = mx_vector)$lt
#' 
life.table <- function(x, Dx = NULL, Ex = NULL, 
                        mx = NULL, qx = NULL, lx0 = 1e+05, ax0 = 0.1){
     nmax     <- length(x)
     n        <- rep(1,nmax)           # width of the intervals
     ax       <- n/2
     if (min(x) == 0) ax[1] <- ax0
     
     mx	   <- if (length(Dx) > 0) { Dx/Ex } 
     else { 
          if (length(mx) > 0) { mx } 
          else {qx/(n - qx*(n - ax)) 
          }   
     }
     if (mx[nmax] < 0.5 | is.na(mx[nmax])) mx[nmax] = mx[nmax - 1]*1.1 
     # In small populations we could have problems 
     # in estimating a reliable mx at last age in the lifetable
     ax[nmax] <- if (mx[nmax] == 0) 0.5 else 1/mx[nmax]
     
     qx       <- if (length(qx) > 0) {qx} 
     else{n*mx / (1 + (n - ax)*mx)} 
     qx[x >= 100 & mx == 0] <- 1
     qx[nmax] <- 1
     
     lx       <- c(1,cumprod(1 - qx))*lx0 
     lx       <- lx[1:nmax]
     dx       <- lx*qx
     Lx       <- n*lx - ax*dx
     Lx[nmax] <- ax[nmax]*dx[nmax]
     Lx[is.na(Lx)] <- 0
     Tx       <- rev(cumsum(rev(Lx)))
     # Tx[nmax] <- max(dx[nmax,], Lx[nmax])
     ex       <- Tx/Lx
     ex[is.na(ex)] <- 0
     ex[nmax] <- if (ex[nmax - 1] == 0) 0 else ax[nmax]
     
     lt <- data.frame(age = x, mx = round(mx,6), qx = round(qx,6), ax = ax, 
                      lx = round(lx), dx = round(dx), Lx = round(Lx), 
                      Tx = round(Tx), ex = round(ex,2))
     lt.exact <- data.frame(age = x, mx = mx, qx = qx, ax = ax,
                            lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex)
     out <- list(lt = lt, lt.exact = lt.exact, process_date = date())
     return(out)
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
     LT <- life.table(ages, mx = fx)
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
#' @return Results
#' @export
#' @examples 
#' library(LinearLink)
#' library(dplyr)
#' 
#' # Select the 1965 - 1990 time interval and fit the Linear-Link model
#' ages    <- 0:100 # available ages in our datasets
#' years   <- 1965:1990 # available years
#' fit_SWE <- HMD_mx$SWE %>% select(mx.1965:mx.1990) %>%
#'      LinearLink(mx = ., ages, years, 'SWEDEN', x_fit = 0)
#' summary(fit_SWE)
#' 
#' # Derive a mortality curve (life table) from a value of
#' # life expectancy at birth in 2014, say 84.05
#' new_e0   <- 84.05
#' pred_SWE <- predict(fit_SWE, new_e0)
#' pred_SWE$lt
LinearLink <- function(mx, mx_ages, mx_years, mx_country = NA, 
                           x_fit = 0, use.smooth = TRUE){
     ptm <- proc.time() # Start the clock!
     
     model_name <- "Linear-Link (2016): ln[m(x)] = b(x)ln[e(x)] + kv(x)"
     # Data preparation
     mx_input <- as.matrix(mx)
     
     # Compute life expectancy
     Life_Tables <- data.frame()
     for (i in 1:ncol(mx)) {
          LT_i   <- life.table(x = mx_ages, mx = mx_input[, i])$lt
          year_i <- mx_years[i]
          LT_i   <- cbind(country = mx_country, year = year_i, LT_i)
          Life_Tables <- rbind(Life_Tables, LT_i)
     }
     colnames(Life_Tables)[3] <- 'age'
     Life_Tables <- Life_Tables[complete.cases(Life_Tables), ] 
     ex  <- Life_Tables[Life_Tables$age == x_fit, ]$ex
     
     # Fit linear portion of model
     x.f <- log(ex)
     y.f <- t(log(mx_input))
     bifit.f <- FUN.bifit(x = x.f, y = y.f)
     
     # Compute residuals and fit SVD portion of model
     yhat1.f <- x.f %*% bifit.f$coef
     dimnames(yhat1.f) <- dimnames(y.f)
     resid1.f <- yhat1.f - y.f
     resid1.f[resid1.f == Inf] <- unique(sort(resid1.f, decreasing = T))[2]
     vx.f  <- svd(resid1.f, 1, 1)$v
     
     # Coefficients -----
     coefs.raw <- round(data.frame(cbind(t(bifit.f$coef), vx.f)), 8)
     coefficients <- coefs.raw
     # I have to smooth the vx's so that we can avoid jumps 
     # in the mortality rates from one age to another. We can do this
     # by using splines. We can allocate 1 degree of freedom for every 5 ages.
     degrees   <- ifelse(use.smooth, round(length(mx_ages)/5), mx_ages)
     df_spline <- ifelse(use.smooth, degrees, 'Smooting not used')
     
     if (use.smooth) {
       coefs.smooth <- coefs.raw*0
       for (j in 1:ncol(coefs.raw)) {
            coefs.smooth[, j] <- smooth.spline(coefs.raw[, j], df = degrees)$y
       }
       coefs.smooth[1, 1] <- coefs.raw[1, 1] # leave infant mortality unsmoothed.
       coefficients <- coefs.smooth
     }
     
     dimnames(coefficients) <- list(mx_ages, c("bx", "vx"))
     #------------------
     
     # Compute fitted values of mx ---- 
     table_ex <- Life_Tables[Life_Tables$age == x_fit, c('year', 'ex')]
     
     LT_optim <- NULL
     k_values <- NULL
     for (i in 1:length(mx_years)) {
          year_i <- table_ex[i, 1]
          ex_target_i <- table_ex[i, 2]
          Optim_out <- FUN.lt_optim(ages = mx_ages, coefs = coefficients, 
                                    ex0 = ex_target_i)
          LT_optim_i <- Optim_out$lt
          LT_optim_i   <- cbind(country = mx_country, 
                                year = year_i, LT_optim_i)
          colnames(LT_optim_i)[3] <- 'age'
          LT_optim <- rbind(LT_optim, LT_optim_i)
          k_optim  <- data.frame(country = mx_country, year = year_i, 
                                 ex = LT_optim_i[LT_optim_i$age == x_fit, ]$ex, 
                                 k = round(Optim_out$k, 6))
          k_values <- rbind(k_values, k_optim)
     }
     table_mx <- LT_optim[, 1:4]
     fitted.values <- reshape(table_mx, direction = 'wide', 
                      idvar = c('country','age'), timevar = 'year')
     fitted.values <- fitted.values[, -(1:2)]
     rownames(fitted.values) <- mx_ages
     residuals <- mx_input - fitted.values
     #-----------------------------------
     
     proc_speed <- round((proc.time() - ptm)[3], 2) # Stop the clock
     cat(paste(mx_country, "Process time:", proc_speed, "seconds!\n") ) # Print speed
     
     out = structure(class = 'LinearLink',
                     list(input_mx = mx_input, input_ages = mx_ages, 
                          input_years = mx_years, input_country = mx_country, 
                          x_fit = x_fit, df_spline = df_spline, 
                coefficients = coefficients, k_values = k_values, 
                fitted.values = fitted.values, residuals = residuals,
                fitted.life.tables = LT_optim, model_name = model_name,
                process_speed = proc_speed, process_date = date()))
     out$call <- match.call()
     return(out)
}









