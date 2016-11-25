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
predict.LinearLink <- function(object, ex_target, ...) {
  coefs <- data.frame(coef(object)$bx, coef(object)$vx, 
                      row.names = object$input$mx_ages)
  pred.values <- FUN.lt_optim(ages = object$input$mx_ages, 
                             coefs = coefs, ex0 = ex_target)
  return(pred.values)
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




