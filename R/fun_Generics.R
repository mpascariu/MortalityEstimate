# S Functions
# Classes and methods

#' @keywords internal
#' @export
summary.LinearLink <- function(object, ...) {
     cat('\nModel:\n')
     cat(object$model_info)
     cat("\n\nCall:\n")
     print(object$call)
     cat('\nDeviance Residuals:\n')
     print(round(summary(as.vector(as.matrix(object$residuals))), 5))
     cat("\nCoefficients:\n")
     coefs <- data.frame(bx = coef(object)$bx, vx = coef(object)$vx,
                         row.names = object$input$mx_ages)
     Hbxvx <- headTail(coefs, digits = 5, hlength = 6, tlength = 6)
     Hk <- headTail(data.frame(. = '.', k = coef(object)$k),
                    digits = 5, hlength = 6, tlength = 6)
     print(data.frame(Hbxvx, Hk))
     cat("\nSpline Smoothing (degrees of freedom): ")
     cat(object$df_spline)
}

#' @keywords internal
#' @export
print.LinearLink <- function(x, ...){
  with(x$input,
       {
        cat('\nLinear-Link Model (2016): \nln[m(x)] = b(x)ln[e(x)] + kv(x)\n')
        cat('\nFitted for life expectancy at age:', theta)
        cat('\nTime interval:', min(mx_years), '-', max(mx_years))
        cat('\nAge-range:', min(mx_ages), '-', max(mx_ages))
        cat('\nCountry:', mx_country, '\n')
        met <- ifelse(method == 'LSE', 'Least Squares (LSE)', 
                      'Poisson Maximum Likelihood (MLE)')
        cat('\nFitting Procedure:', met)
        cat('\nSmoothing:', use.smooth)
  })
}
  
  
#' @keywords internal
#' @export
predict.LinearLink <- function(object, ex_target, 
                               use.vx.rotation = FALSE, ...) {
  # Choose vx coefficients
  vx_o <- coef(object)$vx # original vx coefficients
  vx_r <- rotated_vx(object, e0_t = ex_target) # rotated vx coeffs
  if (use.vx.rotation == TRUE) {vx = vx_r} else {vx = vx_o}
  # Data.frame with all coefficients used in prediction
  coefs <- data.frame(bx = coef(object)$bx, vx = vx, 
                      row.names = object$input$mx_ages)
  # Find the right life table
  pred.values <- FUN.lt_optim(ages = object$input$mx_ages, 
                              coefs = coefs, ex0 = ex_target)
  pred.values$bx <- coefs$bx
  pred.values$vx <- coefs$vx
  return(pred.values)
}

#' @keywords internal
#' 
rotated_vx <- function(object, e0_t, e0_u = 102, 
                       e0_threshold = 80, p_ = 0.5){
  vx = coef(object)$vx
  x  = object$input$mx_ages
  x1 = max(min(x), 15):65
  x2 = 66:max(x)
  
  vx_u     = vx * 0
  vx_young = mean(vx[x1 + 1])
  n_vx     = length(vx[x2 + 1])
  # Derive a logistic shape. The values have to be between 0 and 1, 
  # they will be scaled.
  x_num = seq(-6, 6, length.out = n_vx)
  logit_shape = 1 - exp(x_num)/ (1 + exp(x_num)) 
  vx_old_ = logit_shape * vx_young # scale values
  # This is our new vx
  vx_u[min(x):max(x1 + 1)] = vx_young  # we use only the values between age 15 and 65
  vx_u[x2 + 1] <- vx_old_
  
  # Compute weights
  w_t  <- (e0_t - e0_threshold)/(e0_u - e0_threshold)
  ws_t <- (0.5 * (1 + sin(pi/2 * (2*w_t - 1))) ) ^ p_ 
  # The power to the smooth-weight function, p_, takes values
  # between 0 and 1, which makes the rotation faster at 
  # starting times and slower at ending times. 
  # In this article, p = .5 is taken as the default. (Li et al. 2013)
  
  # Compute rotate_vx
  if ( e0_t < e0_threshold ) {rot_vx = vx }
  if ( e0_threshold <= e0_t & e0_t < e0_u ) {
    rot_vx = (1 - ws_t) * vx + ws_t * vx_u
  }
  if (e0_t >= e0_u) {rot_vx = vx_u}
  return(rot_vx)
}

# ==========================================================================

#' @keywords internal
#' @export
summary.Kannisto <- function(object, ...) {
  cat('Model:\n')
  cat(object$model_name,'\n-----')
  cat("\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  print(headTail(coef(object), digits = 4))
}

#' @keywords internal
#' @export
predict.Kannisto <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) { pred.values <- fitted(object)
  }else{
    x <- newdata
    x_scaled <- x - min(object$x) 
    pars <- coef(object)
    pred.values <- matrix(NA, nrow = length(x), ncol = nrow(pars))
    dimnames(pred.values) <- list(x, rownames(pars))
    fun_ux <- Fun_ux('kannisto')
    for (i in 1:nrow(pars)) {pred.values[,i] = fun_ux(pars[i,], x_scaled)}
  }
  return(pred.values)
}




