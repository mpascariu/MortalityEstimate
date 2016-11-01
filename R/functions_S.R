# S Functions
# Classes and methods

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
LinearLink <- function(mx, mx_ages, mx_years, mx_country = NA, x_fit = 0, 
                       use.smooth = TRUE) UseMethod("LinearLink")

#' @keywords internal
#' @export
LinearLink.default <- function(mx, mx_ages, mx_years, mx_country = NA, 
                               x_fit = 0, use.smooth = TRUE){
     mdl <- FUN.linearlink(mx, mx_ages, mx_years, mx_country, x_fit, use.smooth)
     mdl$call   <- match.call()
     class(mdl) <- "LinearLink"
     mdl
}

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






