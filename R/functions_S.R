# S Functions
# Classes and methods

#' @keywords internal
#' @export
summary.LinearLink <- function(object, ...) {
     cat('Model:\n')
     cat(object$model_name,'\n-----')
     cat("\nCall:\n")
     print(object$call)
     cat("\nCoefficients:\n")
     print(headTail(coef(object), digits = 5))
     cat("\nk-values:\n")
     print(headTail(object$k_values, digits = 6))
     cat("\n\nSpline Smoothing (degrees of freedom): ")
     cat(object$df_spline)
}


#' @keywords internal
#' @export
predict.LinearLink <- function(object, ex_target, ...) {
     pred.values <- FUN.lt_optim(ages = object$input_ages, 
                                 coefs = coef(object), ex0 = ex_target)
     pred.values
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
  if(is.null(newdata)){pred.values <- fitted(object)
  }else{
    x <- newdata
    x_scaled <- x - min(object$x) 
    pars <- coef(object)
    pred.values <- matrix(NA, nrow = length(x), ncol = nrow(pars))
    dimnames(pred.values) <- list(x, rownames(pars))
    fun_ux <- Fun_ux('kannisto')
    for(i in 1:nrow(pars)) pred.values[,i] = fun_ux(pars[i,], x_scaled)
  }
  pred.values
}




