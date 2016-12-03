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
  

#' Predict function for a LinearLink object
#' 
#' @param object An object of class 'LinearLink'
#' @param use.vx.rotation Logical argument. If \code{TRUE} the adjustment method
#' described in Li et al. (2013) paper is applied to the vx coefficients before 
#' estimated the life table. If \code{FALSE} the fitted vx coefficients are used 
#' in the estimations of the life table.
#' @param ex_target A value of life expectancy for which we want to estimate 
#' the mortality curve
#' @inheritParams rotated_vx 
#' @param ... Other arguments 
# #' @keywords internal
#' @source Li et al. (2013) 
#' @export
predict.LinearLink <- function(object, ex_target, 
                               use.vx.rotation = FALSE, ...) {
  # Choose vx coefficients
  if (use.vx.rotation == TRUE) { 
    vx = rotated_vx(object, ex_target = ex_target, ...) } else { 
      vx = coef(object)$vx }
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


#' Compute rotated vx coefficients
#' 
#' This functions computes the rotation of the vx coefficients using the method
#' presented in Li et al. (2013) paper. 
#' @param object An object of class 'LinearLink'
#' @param ex_target A value of life expectancy for which we want to derive the rotated vx
#' @param e0_u Ultimate value of life expectancy. At this point the rotation 
#' process reaches its maximum efficiency. Here. e0_u = 80 is taken as the default.  
#' @param e0_threshold Level of life expectancy where the rotation should begin.
#' If rotated_vx is computed for ex_target <= e0_threshold then no diffrence
#' will be observed.  
#' @param p_ The power to the smooth-weight function, p_, takes values 
#' between 0 and 1, which makes the rotation faster at starting times and 
#' slower at ending times. Here, p = .5 is taken as the default
#' @param positive.vx.only Logical argument. If negative values of vx are 
#' estimated we can shift up the entire vx curve so that the minimum value 
#' is equal to zero.
#' @source Li et al. (2013) 
#' @keywords internal
#' 
rotated_vx <- function(object, ex_target, e0_u = 102, 
                       e0_threshold = 80, p_ = 0.5,
                       positive.vx.only = FALSE){
  vx = coef(object)$vx
  x  = object$input$mx_ages
  x1 = max(min(x), 15):65 # young ages
  x2 = 66:max(x) # old ages
  
  vx_u     = vx * 0
  vx_young_ages = mean(vx[x1 + 1]) # select vx corresponding to young ages. 
  # Only the values between age 15 and 65 are being used.
  
  n_vx     = length(vx[x2]) #count age groups in x2
  # Derive a logistic shape. The values have to be between 0 and 1, 
  # they will be scaled.
  x_num = seq(-6, 6, length.out = n_vx)
  logit_shape = 1 - exp(x_num)/ (1 + exp(x_num)) 
  vx_old_ages = logit_shape * vx_young_ages # scale values
  # This is our vx ultimate (vx_u)
  vx_u[min(x):max(x1 + 1)] <- vx_young_ages
  vx_u[x2 + 1] <- vx_old_ages
  # Compute weights
  w_t  = (ex_target - e0_threshold)/(e0_u - e0_threshold)
  ws_t = (0.5 * (1 + sin(pi/2 * (2*w_t - 1))) ) ^ p_ 
  # Compute rotate_vx
  if ( ex_target < e0_threshold ) { rot_vx = vx }
  if ( e0_threshold <= ex_target & ex_target < e0_u ) {
    rot_vx = (1 - ws_t) * vx + ws_t * vx_u
  }
  if (ex_target >= e0_u) { rot_vx = vx_u }
  # Shift up vx curve (suggested by Ugo Basellini)
  if (positive.vx.only == TRUE & min(rot_vx) < 0) { 
    rot_vx = rot_vx + abs(min(rot_vx)) 
    }
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




