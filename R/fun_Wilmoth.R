
#' Fit Log-Quadratic model
#' 
#' Estimating the log-quadratic model using the bi-weight procedure as 
#' described in the Appendix of Wilmoth et.al.(2012).
#' @param data input data. A collection of life tables.
#' @param sex Specify the sex of the population. Options: \code{"female", "male" or "total"}
#' @param ages Numerical vector containing ages covered in the life tables.
#' If abridge life tables are supplied in \code{data} then provide the lower 
#' bounds of the age intervals. Default: \code{c(0, 1, seq(5, 110, by = 5))} 
#' (corresponding to the HMD format).
#' @return The output is of class \code{wilmoth} with the components:
#' @return \item{input}{ a list with input objects provided in the function}
#' @return \item{coefficients}{ estimated coefficient}
#' @return \item{fitted.values}{ fitted values}
#' @return \item{model.life.table}{ a model life table}
#' @return \item{model_info}{ Model details (equqtion). The relationship between 
#' the death rate at age x, mx, and the probability of dying between birth and 
#' age 5.}
#' @references John Wilmoth, Sarah Zureick, Vladimir Canudas-Romo, Mie Inoue & 
#' Cheryl Sawyer (2012): A flexible two-dimensional mortality model for use in 
#' indirect estimation, Population Studies: A Journal of Demography, 66:1, 1-28
#' \url{http://dx.doi.org/10.1080/00324728.2011.611411}
#' @seealso \code{\link{wilmothLT}}
#' @examples 
#' 
#' # DATA
#' sex = "female"
#' HMD719f <- HMD719[HMD719$sex == sex, ]
#' 
#' # Fit model
#' W <- wilmoth(data = HMD719f, sex)
#' 
#' # Build life tables with various choices of 2 input parameters ---
#' # (case 1) Using 5q0 and k
#' WLT_1 <- wilmothLT(W, q0_5 = 0.05, k = 0.1)
#' WLT_1
#' @export  
#'  
wilmoth <- function(data, sex, ages = c(0,1,seq(5, 110, by = 5))) {
  
  if (!(sex %in% c("female", "male", "total"))) stop("sex must be: 'female', 'male', or 'total'", call. = FALSE)
  
  info    <- "log m[x] = a[x] + b[x] h + c[x] h^2 + kv[x]"
  input    <- c(as.list(environment()))
  na       <- length(ages)
  gr_names <- paste0("[", ages,",", c(ages[-1], "+"), ")")
  data$gr_names <- gr_names
  
  # 5q0 (vector)
  q0_5 <- 1 - data[data$gr_names == "[5,10)", "lx"] / data[data$gr_names == "[0,1)", "lx"]
  # mx (matrix)
  mx_ <- array(data$mx, dim = c(na, nrow(data)/na))
  
  # Fit log-quadratic portion of model
  x.f     <- cbind(log(q0_5), log(q0_5)^2)
  y.f     <- t(log(mx_))
  bifit.f <- bifit(x.f, y.f, c = 6)
  coef    <- data.frame(t(bifit.f$coef))
  dimnames(coef) <- list(gr_names, c("ax", "bx","cx"))
  # coef[2, ] <- NA
  
  # Compute residuals and fit SVD portion of model
  yhat1.f  <- cbind(1, x.f) %*% t(coef)
  dimnames(yhat1.f) <- dimnames(y.f)
  resid1.f <- yhat1.f - y.f
  svd.f    <- svd(resid1.f[, -c(1,2)], 1, 1)
  # Set values of vx = 0 for ages 0, 1-4 and above 90
  vx <- c(0, NA, svd.f$v) 
  vx[20:24] <- 0
  coef$vx   <- vx
  
  H <- median(q0_5)
  fit_mx <- mxhat.logquad(coef, x = ages, sex, q0_5 = H, k = 0)
  fit_lt <- lthat.logquad(coef, x = ages, sex, q0_5 = H, k = 0)$lt
  
  out <- list(input = input, coefficients = coef, 
              fitted.values = fit_mx, model.life.table = fit_lt, 
              model.info = info)
  out <- structure(class = 'wilmoth', out)
  out$call <- match.call()
  return(out)
}


#' Estimate wilmoth life table
#' 
#' Construct a life table based on the log-quadratic (wilmoth) estimates
#' with various choices of 2 input parameters:  
#' \code{5q0, k, e0, 45q15, 35q15, 1q0}.
#' @param object An \code{\link{wilmoth}} object.
#' @param q0_5 5q0. The probability that a newborn will die during the subsequent 5 years
#' @param k k-parameter in the log-quadratic model
#' @param e0 life expectancy at birth
#' @param q15_45 45q15. The probability that a life aged 15 will die during the subsequent 45 years
#' @param q15_35 35q15. The probability that a life aged 15 will die during the subsequent 35 years
#' @param q0_1 1q0. The probability that a life aged 0 will die during the following year
#' @param tol tolerance level for convergence. The tolerance level, is relevant 
#' for case 12 and 13 (QQ and e0 are known)
#' @param maxit maximum number of iterations allowed. Default: 100.
#' @param ... additional arguments affecting the predictions produced.
#' @return a life table matching given inputs, and associated values of 
#' \code{5q0, k, e0, 45q15, 35q15,} and \code{1q0}.
#' @seealso \code{\link{wilmoth}}
#' @examples 
#' # checck wilmoth examples
#' @export
#' 
wilmothLT <- function(object, q0_5 = NULL, k = NULL, e0 = NULL, 
                            q15_45 = NULL, q15_35 = NULL, q0_1 = NULL, 
                            tol = 1e-9, maxit = 200, ...) {
  
  my_case <- find.my.case(q0_5, k, e0, q15_45, q15_35, q0_1)
  cf  <- coef(object)
  sex <- object$input$sex
  x   <- object$input$ages
  
  # Cases 1-4:  5q0 is known, plus k, e0, 45q15 or 45q15
  if (my_case == "C1") {
    tmp <- lthat.logquad(cf, x, sex, q0_5, k)
  }
  if (my_case %in% c("C2", "C3", "C4")) {
    if (my_case == "C2") fun.k <- function(k) {
      lthat.logquad(cf, x, sex, q0_5, k)$lt$ex[1] - e0
    }
    if (my_case == "C3") fun.k <- function(k) { 
      lt = lthat.logquad(cf, x, sex, q0_5, k)$lt
      return(1 - lt["[60,65)", "lx"] / lt["[15,20)", "lx"] - q15_45)
    }
    if (my_case == "C4") fun.k <- function(k) { 
      lt = lthat.logquad(cf, x, sex, q0_5, k)$lt
      return(1 - lt["[50,55)", "lx"] / lt["[15,20)", "lx"] - q15_35)
    }
    
    root <- uniroot(f = fun.k, interval = c(-10,10))$root
    tmp  <- lthat.logquad(cf, x, sex, q0_5, k = root) 
  }
  
  # Cases 5-8:  
  # 1q0 is known, plus k, e0, 45q15 or 35q15;
  # after finding 5q0 (assume k=0, but it doesn't matter), these become Cases 1-4
  if (my_case %in% c("C5","C6","C7","C8") ) {
    fun.q0_5  <- function(q0_5) lthat.logquad(cf, x, sex, q0_5, k = 0)$lt$qx[1] - q0_1
    root <- uniroot(f = fun.q0_5, interval = c(1e-5, 0.8))$root
  }
  if (my_case == "C5") tmp <- wilmothLT(object, q0_5 = root, k = k, tol = tol, maxit = maxit, ...)
  if (my_case == "C6") tmp <- wilmothLT(object, q0_5 = root, e0 = e0, tol = tol, maxit = maxit, ...)
  if (my_case == "C7") tmp <- wilmothLT(object, q0_5 = root, q15_45 = q15_45, tol = tol, maxit = maxit, ...)
  if (my_case == "C8") tmp <- wilmothLT(object, q0_5 = root, q15_35 = q15_35, tol = tol, maxit = maxit, ...)
  
  # Cases 9-11:  k is known, plus e0, 45q15 or 35q15; 
  # must find 5q0
  if (my_case %in% c("C9", "C10", "C11")) {
    if (my_case == "C9") fun.q0_5 = function(q0_5) { 
      lthat.logquad(cf, x, sex, q0_5, k)$lt$ex[1] - e0 
    }
    if (my_case == "C10") fun.q0_5 = function(q0_5) { 
      lt <- lthat.logquad(cf, x, sex, q0_5, k)$lt
      return(1 - lt["[60,65)", "lx"] / lt["[15,20)", "lx"] - q15_45)
    }
    if (my_case == "C11") fun.q0_5 <- function(q0_5) {
      lt <- lthat.logquad(cf, x, sex, q0_5, k)$lt
      return(1 - lt["[50,55)", "lx"] / lt["[15,20)", "lx"] - q15_35)
    }
    
    root <- uniroot(fun.q0_5, c(1e-4,0.8))$root
    tmp  <- lthat.logquad(cf, x, sex, q0_5 = root, k)
  }
  
  # Case 12 and 13: QQ and e0 are known; must find both 5q0 and k
  if (my_case %in% c("C12", "C13")) {
    k <- q0_5 <- 0
    iter <- crit <- 1
    
    while (crit > tol & iter <= maxit) {
      k.old  <- k
      q0_5.old <- q0_5
      q0_5 <- predict(object, k = k, e0 = e0, tol = tol, maxit = maxit, ...)$q0_5 # Get 5q0 from e0 assuming k
      # Get k from 45q15 or 35q15 asuming 5q0
      if (my_case == "C12") tmp = wilmothLT(object, q0_5 = q0_5, q15_45 = q15_45, tol = tol, maxit = maxit, ...)
      if (my_case == "C13") tmp = wilmothLT(object, q0_5 = q0_5, q15_35 = q15_35, tol = tol, maxit = maxit, ...)
      k  <- tmp$k
      crit = sum(abs(c(k, q0_5) - c(k.old, q0_5.old)))
      iter = iter + 1
    }
    if (iter > maxit) warning("number of iterations reached maximum without convergence", call. = F)
  }
  
  # Extract lt, lt.exact, 5q0, and k from tmp
  lt   <- tmp$lt
  q0_5 <- tmp$q0_5
  k    <- tmp$k
  # Extract e0, 45q15, 35q15, and 1q0 from final (exact) life table
  e0     <- lt$ex[1]
  q0_1   <- lt$qx[1]
  q15_45 <- 1 - lt["[60,65)", "lx"] / lt["[15,20)", "lx"]
  q15_35 <- 1 - lt["[50,55)", "lx"] / lt["[15,20)", "lx"]
  
  # mx1x1 <- mxfun(coefs, sex, q0_5, k, x = 0.5:110.5) 
  # names(mx1x1) <- 0:110
  
  # Return life table plus values of the 6 possible inputs
  out = list(lt = lt, indices = data.frame(k, q0_1, q0_5, q15_35, q15_45, e0))
  out = structure(class = "wilmothLT", out)
  return(out)
}

#' Compute WLS using Tukey's biweight function
#' 
#' @param x a matrix whose rows correspond to cases and whose columns correspond to variables.
#' x matrix should not include a column of ones
#' @param y the responses, possibly a matrix if you want to fit multiple left hand sides.
#' y can be matrix if there are multiple left-hand sides
#' @param c the tuning constant, typically around 6 to 9
#' @param tol the tolerance for convergence of the biweight fit
#' @param intercept logical. Whether or not an intercept term should be used.
#' @param tolerance the tolerance to be used in the matrix decomposition;
#' @param yname names to be used for the response variables.
#' @keywords internal
#' 
bifit <- function(x, y, c = 6, tol = 1e-6, 
                  intercept = TRUE, tolerance = 1e-07, yname = NULL) {
	
	if (is.vector(y)) {
		coef.old <- iter <- 0
		coef <- 1 
		wt <- NULL
		while (max(abs(coef - coef.old)) > tol) {
			iter <- iter + 1
			z    <- lsfit(x, y, wt, intercept = intercept, tolerance = tolerance, yname = yname)
			S    <- median(abs(z$resid))
			u    <- z$resid / (c*S)
			wt   <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
			coef.old <- coef
			coef <- z$coef 
		}
		resid <- z$resid
		names(resid) <- names(wt) <- names(u) <- names(y)
		names(coef) <- paste0("b", 0:(ncol(cbind(1, x)) - 1))
	}

	if (is.matrix(y)) {
		resid <- wt <- u <- coef <- NULL
		for (j in 1:ncol(y)) {
			z <- bifit(x, y[,j], c = c, tol = tol, intercept = intercept, 
			           tolerance = tolerance, yname = yname)
			resid <- cbind(resid, z$resid)
			wt <- cbind(wt, z$wt)
			u <- cbind(u, z$u)
			coef <- cbind(coef, z$coef) }
		dimnames(resid) <- dimnames(wt) <- dimnames(u) <- dimnames(y)
		dimnames(coef) <- list(paste0("b", 0:(ncol(cbind(1, x)) - 1)), dimnames(y)[[2]])
	}
  out <- list(coef = coef, residuals = resid, intercept = intercept, wt = wt, u = u)
  return(out)
}


#' Estimate mx using the log-quadratic model
#' @inheritParams lthat.logquad
#' @keywords internal
#' 
mxhat.logquad <- function(coefs, x, sex, q0_5, k) {
  h  <- log(q0_5)
	mx <- with(as.list(coefs), exp(ax + bx*h + cx*h^2 + vx*k))
	qx <- mx_qx(x, mx, out = "qx")
	
	# Force 4q1 (and thus 4m1) to be consistent with 1q0 and 5q0
	a4    <- coale.demeny.ax(x, mx, qx, sex)[2]
	q0_1  <- mx[1] / (1 + (1 - a4)*mx[1])
	q1_4  <- 1 - (1 - q0_5)/(1 - q0_1)
	mx[2] <- q1_4 / (4 - (4 - a4)*q1_4)
	
	names(mx) <- rownames(coefs)
	return(mx) 
}


#' Estimated life table using the log-quadratic model
#' 
#' @param x numeric vector of ages.
#' @param mx numeric vector of death dates
#' @param sex sex. Possible values: \code{"female", "male"}, or \code{"total"}.
#' @param q0_5 5q0
#' @param lx0 radix. Default: 100 000.
#' @keywords internal
#' 
lthat.logquad <- function(coefs, x, sex, q0_5, k, lx0 = 10^5) {
  mx <- mxhat.logquad(coefs, x, sex, q0_5, k)
  nx <- c(diff(x), Inf)
  qx <- mx_qx(x, mx, out = "qx")
  ax <- coale.demeny.ax(x, mx, qx, sex) # Below age 5
  
  # 5q0 is given. 
  # 4q1 is derived from 1q0 and 5q0 => 4m1 is re-computed 
  if (!is.null(q0_5)) {
    qx[2] <- 1 - (1 - q0_5)/(1 - qx[1])
    mx[2] <- qx[2] / (nx[2] - (nx[2] - ax[2])*qx[2] ) 
  }
  
  LT  <- lt.core(x, mx, qx, ax, lx0)$lt.exact
  out <- list(lt = LT, q0_5 = q0_5, k = k)
  return(out)
}


#' Function that determintes the case/problem we have to solve
#' It also performes some checks
#' @keywords internal
#' 
find.my.case <- function(q0_5, k, e0, q15_45, q15_35, q0_1) {
  # Test that at least of one of 1q0 and 5q0 is null
  my_case <- check.is.null(q0_5, k, e0, q15_45, q15_35, q0_1)
  
  if (sum(my_case[c(1, 6)]) == 0) stop("cannot have both 'q0_1' and 'q0_5' as inputs", call. = FALSE)
  # Test that at least of one of 45q15 and 35q15 is null
  if (sum(my_case[c(4, 5)]) == 0) stop("cannot have both 'q15_45' and 'q15_35' as inputs", call. = FALSE)
  # Test that exactly two inputs are non-null
  if (sum(my_case) != 4) stop("must have exactly two inputs", call. = FALSE)
  # There are 13 cases:  "6 choose 2" = 15, but we disallow two cases (1q0 and 5q0, or 45q15 and 35q15)
  if (all(my_case == c(F,F,T,T,T,T))) case = "C1" 
  if (all(my_case == c(F,T,F,T,T,T))) case = "C2" 
  if (all(my_case == c(F,T,T,F,T,T))) case = "C3" 
  if (all(my_case == c(F,T,T,T,F,T))) case = "C4"
  if (all(my_case == c(T,F,T,T,T,F))) case = "C5" 
  if (all(my_case == c(T,T,F,T,T,F))) case = "C6"   
  if (all(my_case == c(T,T,T,F,T,F))) case = "C7" 
  if (all(my_case == c(T,T,T,T,F,F))) case = "C8" 
  if (all(my_case == c(T,F,F,T,T,T))) case = "C9" 
  if (all(my_case == c(T,F,T,F,T,T))) case = "C10" 
  if (all(my_case == c(T,F,T,T,F,T))) case = "C11" 
  if (all(my_case == c(T,T,F,F,T,T))) case = "C12" 
  if (all(my_case == c(T,T,F,T,F,T))) case = "C13" 
  return(case)
}


#' @keywords internal
#' 
check.is.null <- function(...) {
  ls = list(...)
  unlist(lapply(ls, is.null))
}

# ----------------------------------------------
# S3 functions

#' @keywords internal
#' @export
print.wilmoth <- function(x, ...){
  cat("\nThe Log-Quadratic model:\n")
  cat(x$model.info, "\n")
}

#' @keywords internal
#' @export
summary.wilmoth <- function(object, ...) print.wilmoth(x = object, ...)