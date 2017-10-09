
#' Fit log-quadratic model using bi-weight method
#' 
#' @param data input data. A collection of life tables.
#' @param sex Specify the sex of the population. Options: \code{"Female", "Male" or "Total"}
#' @param ages Numerical vector containing ages covered in the life tables.
#' If abridge life tables are supplied in \code{data} then provide the lower 
#' bounds of the age intervals. Default: \code{c(0, 1, seq(5, 110, by = 5))} 
#' (corresponding to the HMD format).
#' @return An \code{wilmoth} object
#' @export  
#'  
wilmoth <- function(data, sex, ages = c(0,1,seq(5, 110, by = 5))) {
  
  if (!(sex %in% c("Female", "Male", "Total"))) stop("sex must be: 'Female', 'Male', or 'Total'", call. = FALSE)
  
  input <- c(as.list(environment()))
  na        <- length(ages)
  ageG      <- paste0("[", ages,",", c(ages[-1], "+"), ")")
  ageG[na]  <- paste0("[", ages[na], "+]")
  data$ageG <- ageG
  
  # 5q0 (vector)
  q0_5 <- 1 - data[data$ageG == "[5,10)", "lx"] / data[data$ageG == "[0,1)", "lx"]
  # mx (matrix)
  mx_ <- array(data$mx, dim = c(na, nrow(data)/na))
  
  # Fit log-quadratic portion of model
  x.f     <- cbind(log(q0_5), log(q0_5)^2)
  y.f     <- t(log(mx_))
  bifit.f <- bifit(x.f, y.f, c = 6)
  coef    <- data.frame(t(bifit.f$coef))
  dimnames(coef) <- list(ageG, c("ax", "bx","cx"))
  # coef[2, ] <- NA
  
  # Compute residuals and fit SVD portion of model
  yhat1.f  <- cbind(1, x.f) %*% t(coef)
  dimnames(yhat1.f) <- dimnames(y.f)
  resid1.f <- yhat1.f - y.f
  svd.f    <- svd(resid1.f[, -c(1,2)], 1, 1)
  # Set values of vx = 0 for ages 0, 1-4 and above 90
  vx <- c(0, NA, svd.f$v) 
  vx[20:24] <- 0
  coef$vx <- vx
  
  fit_mx <- mxhat.logquad(coef, sex, Q5 = median(q0_5))
  fit_lt <- lt.from.mx(mx = fit_mx, ages, sex, Q5 = median(q0_5))
  
  out <- list(input = input, coefficients = coef, 
              fitted.values = fit_mx, model.life.table = fit_lt,
              ages_groups = ageG)
  out <- structure(class = 'wilmoth', out)
  out$call <- match.call()
  return(out)
}


#' Predict function for wilmoth object
#' 
#' @param object An \code{\link{wilmoth}} object.
#' @param Q5 5q0
#' @param k k 
#' @param e0 life expectancy at birth level
#' @param QQa 45q15
#' @param QQb 35q15
#' @param Q1 1q0
#' @param tol tolerance level for convergence. The tolerance level, is relevant 
#' for case 12 and 13 (QQ and e0 are known)
#' @param maxit maximum number of iterations allowed. Default: 100.
#' @param ... additional arguments affecting the predictions produced.
#' @return life table(s) matching given inputs, and associated values of 
#' Q5, k, e0, QQa, QQb, and Q1
#' @export
#' 
predict.wilmoth <- function(object, Q5 = NULL, k = NULL, e0 = NULL, 
                            QQa = NULL, QQb = NULL, Q1 = NULL, 
                            tol = 1e-9, maxit = 200, ...) {
  
  my_case <- find.my.case(Q5, k, e0, QQa, QQb, Q1)
  cf   <- coef(object)
  sex  <- object$input$sex
  ages <- object$input$ages
  
  # Cases 1-4:  Q5 is known, plus k, e0, QQa or QQb
  if (my_case == "C1") {
    tmp <- lthat.logquad(cf, ages, sex, Q5, k)
  }
  if (my_case %in% c("C2", "C3", "C4")) {
    if (my_case == "C2") fun.k <- function(k) {
      lthat.logquad(cf, ages, sex, Q5, k)$lt$ex[1] - e0
    }
    if (my_case == "C3") fun.k <- function(k) { 
      lt = lthat.logquad(cf, ages, sex, Q5, k)$lt
      return(1 - lt["[60,65)", "lx"] / lt["[15,20)", "lx"] - QQa)
    }
    if (my_case == "C4") fun.k <- function(k) { 
      lt = lthat.logquad(cf, ages, sex, Q5, k)$lt
      return(1 - lt["[50,55)", "lx"] / lt["[15,20)", "lx"] - QQb)
    }
    
    root <- uniroot(f = fun.k, interval = c(-10,10))$root
    tmp  <- lthat.logquad(cf, ages, sex, Q5, k = root) 
  }
  
  # Cases 5-8:  Q1 is known, plus k, e0, QQa or QQb;
  #             after finding Q5 (assume k=0, but it doesn't matter), these become Cases 1-4
  if (my_case %in% c("C5","C6","C7","C8") ) {
    fun.Q5  <- function(Q5) lthat.logquad(cf, ages, sex, Q5, k = 0)$lt$qx[1] - Q1
    root <- uniroot(f = fun.Q5, interval = c(1e-5, 0.8))$root
  }
  if (my_case == "C5") tmp <- predict(object, Q5 = root, k = k, tol = tol, maxit = maxit, ...)
  if (my_case == "C6") tmp <- predict(object, Q5 = root, e0 = e0, tol = tol, maxit = maxit, ...)
  if (my_case == "C7") tmp <- predict(object, Q5 = root, QQa = QQa, tol = tol, maxit = maxit, ...)
  if (my_case == "C8") tmp <- predict(object, Q5 = root, QQb = QQb, tol = tol, maxit = maxit, ...)
  
  # Cases 9-11:  k is known, plus e0, QQa or QQb; 
  # must find Q5
  if (my_case %in% c("C9", "C10", "C11")) {
    if (my_case == "C9") fun.Q5 = function(Q5) { 
      lthat.logquad(cf, ages, sex, Q5, k)$lt$ex[1] - e0 
    }
    if (my_case == "C10") fun.Q5 = function(Q5) { 
      lt <- lthat.logquad(cf, ages, sex, Q5, k)$lt
      return(1 - lt["[60,65)", "lx"] / lt["[15,20)", "lx"] - QQa)
    }
    if (my_case == "C11") fun.Q5 <- function(Q5) {
      lt <- lthat.logquad(cf, ages, sex, Q5, k)$lt
      return(1 - lt["[50,55)", "lx"] / lt["[15,20)", "lx"] - QQb)
    }
    
    root <- uniroot(fun.Q5, c(1e-4,0.8))$root
    tmp  <- lthat.logquad(cf, ages, sex, Q5 = root, k)
  }
  
  # Case 12 and 13: QQ and e0 are known; must find both Q5 and k
  if (my_case %in% c("C12", "C13")) {
    k <- Q5 <- 0
    iter <- crit <- 1
    
    while (crit > tol & iter <= maxit) {
      k.old  <- k
      Q5.old <- Q5
      Q5 <- predict(object, k = k, e0 = e0, tol = tol, maxit = maxit, ...)$Q5 # Get Q5 from e0 assuming k
      # Get k from QQa or QQb asuming Q5
      if (my_case == "C12") tmp = predict(object, Q5 = Q5, QQa = QQa, tol = tol, maxit = maxit, ...)
      if (my_case == "C13") tmp = predict(object, Q5 = Q5, QQb = QQb, tol = tol, maxit = maxit, ...)
      k  <- tmp$k
      crit = sum(abs(c(k, Q5) - c(k.old, Q5.old)))
      iter = iter + 1
    }
    if (iter > maxit) warning("number of iterations reached maximum without convergence", call. = F)
  }
  
  # Extract lt, lt.exact, Q5, and k from tmp
  lt  <- tmp$lt
  Q5  <- tmp$Q5
  k   <- tmp$k
  # Extract e0, QQa, QQb, and Q1 from final (exact) life table
  e0  <- lt$ex[1]
  Q1  <- lt$qx[1]
  QQa <- 1 - lt["[60,65)", "lx"] / lt["[15,20)", "lx"]
  QQb <- 1 - lt["[50,55)", "lx"] / lt["[15,20)", "lx"]
  
  # mx1x1 <- mxfun(coefs, sex, Q5, k, x = 0.5:110.5) 
  # names(mx1x1) <- 0:110
  
  # Return life table plus values of the 6 possible inputs
  out <- list(lt = lt, Q5 = Q5, k = k, 
              e0 = e0, QQa = QQa, QQb = QQb, Q1 = Q1)#, mx1x1 = mx1x1)
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
mxhat.logquad <- function(coefs, sex, Q5, k = 0) {
	mx <- with(as.list(coefs), exp(ax + bx*log(Q5) + cx*log(Q5)^2 + vx*k))
	
	# Force 4q1 (and thus 4m1) to be consistent with 1q0 and 5q0
	a4 <- coale.demeny.ax(mx[1], sex)[2]
	Q1 <- mx[1] / (1 + (1 - a4)*mx[1] )
	Q4 <- 1 - (1 - Q5)/(1 - Q1)
	mx[2] <- Q4 / (4 - (4 - a4)*Q4 )
  	
	names(mx) <- rownames(coefs)
	return(mx) 
}


#' Estimated life table using the log-quadratic model
#' @keywords internal
#' 
lthat.logquad <- function(coefs, ages, sex, Q5, k) {
	mxhat  <- mxhat.logquad(coefs, sex, Q5, k)
	lt     <- lt.from.mx(mx = mxhat, ages, sex, Q5)
	out    <- list(lt = lt, Q5 = Q5, k = k)
	return(out) 
}


#' Find ax indicator using the Coale-Demeny coefficients
#' 
#' ax - the point in the age internal where 50% of the deaths have already occurred
#' @keywords internal
#' 
coale.demeny.ax <- function(m0, sex) {
	if (sum(m0 <= 0) > 0) stop("m0 must be greater than 0", call. = F)
  
  a0M   <- ifelse(m0 >= 0.107, 0.330, 0.045 + 2.684*m0)
  a1_4M <- ifelse(m0 >= 0.107, 0.330, 1.651 - 2.816*m0)
  a0F   <- ifelse(m0 >= 0.107, 0.350, 0.053 + 2.800*m0)
  a1_4F <- ifelse(m0 >= 0.107, 0.350, 1.522 - 1.518*m0)
  a0T   <- (a0M + a0F)/2
  a1_4T <- (a1_4M + a1_4F)/2
	  
  out <- matrix(c(a0F, a1_4F, a0M, a1_4M, a0T, a1_4T), ncol = 2, byrow = T)
  dimnames(out) <- list(c("Female", "Male", "Total"), c("a0", "a1_4"))
	return(out[sex, ]) 
}

#' Construct a life table from a vector of death rates
#' 
#' @param mx numeric vector of death dates
#' @param ages numeric vector of ages.
#' @param sex sex. Possible values: \code{"female", "male"}, or \code{"total"}.
#' @param Q5 Q5
#' @param lx0 radix. Default: 100 000.
#' @return 
#' @export
#' 
lt.from.mx <- function(mx, ages, sex, Q5 = NULL, lx0 = 10^5) {
  d_names <- list(Age_Group = names(mx), c("mx","qx","ax","lx","dx","Lx","Tx","ex"))
  
  N  <- length(mx)
  nx <- c(diff(ages), Inf)
  qx <- 1 - exp(-nx*mx)
  ax <- nx + 1/mx - nx/qx 
  # Below age 5
  ax[1:2] <- coale.demeny.ax(mx[1], sex)
  qx[1] <- nx[1]*mx[1] / (1 + (nx[1] - ax[1])*mx[1])
  # If Q5 is given, 4q1 is derived from 1q0 and 5q0 and 4m1 is re-computed; 
  # otherwise, 4q1 is derived directly from 4m1
  if (is.null(Q5)) {
    qx[2] <- nx[2]*mx[2] / (1 + (nx[2] - ax[2])*mx[2])
  } else {
    qx[2] <- 1 - (1 - Q5)/(1 - qx[1])
    mx[2] <- qx[2] / (nx[2] - (nx[2] - ax[2])*qx[2] ) 
  }
  
  # Open age interval
  qx[N] <- 1
  ax[N] <- 1/mx[N]
  # Make life table
  lx <- lx0 * c(1, cumprod(1 - qx)[1:(N - 1)])
  dx <- qx * lx
  Lx <- nx*lx - (nx - ax)*dx 
  Lx[N] <- dx[N]/mx[N]
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx/lx
  lt <- data.frame(mx, qx, ax, lx, dx, Lx, Tx, ex)
  dimnames(lt) <- d_names
  return(lt)
}


#' Function that determintes the case/problem we have to solve
#' It also performes some checks
#' @keywords internal
#' 
find.my.case <- function(Q5, k, e0, QQa, QQb, Q1) {
  # Test that at least of one of Q1 and Q5 is null
  my_case <- check.is.null(Q5, k, e0, QQa, QQb, Q1)
  
  if (sum(my_case[c(1, 6)]) == 0) stop("cannot have both Q1 and Q5 as inputs", call. = FALSE)
  # Test that at least of one of QQa and QQb is null
  if (sum(my_case[c(4, 5)]) == 0) stop("cannot have both QQa = 45q15 and QQb = 35q15 as inputs", call. = FALSE)
  # Test that exactly two inputs are non-null
  if (sum(my_case) != 4) stop("must have exactly two inputs", call. = FALSE)
  # There are 13 cases:  "6 choose 2" = 15, but we disallow two cases (Q1 and Q5, or QQa and QQb)
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