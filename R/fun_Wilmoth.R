# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# Last update: Fri Apr 30 14:07:41 2021
# --------------------------------------------------- #

#' Fit the Log-Quadratic Model
#' 
#' Estimate the Log-Quadratic model. The implemented estimation using the 
#' bi-weight procedure is described in the Appendix of 
#' \insertCite{wilmoth2012;textual}{MortalityEstimate}
#' @param x Numerical vector containing ages corresponding to the input data 
#' (in 'mx' or 'LT').
#' @param mx Input data. A data.frame / matrix with death rates.
#' @param LT Input data. A collection of life tables. If \code{mx} is provided 
#' \code{LT} is not necessary and vice versa.
#' @param verbose Logical. Choose whether to display a progress bar during the 
#' fitting process. Default: TRUE.
#' @param control List with additional parameters. 
#' See \code{\link{wilmoth_control}}.
#' @return The output is of class \code{wilmoth} with the components:
#'  \item{input}{A list with input objects provided in the function;}
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given 
#'  arguments;}
#'  \item{coefficients}{Estimated coefficients;}
#'  \item{k}{The estimated k parameter;}
#'  \item{fitted.values}{Fitted values of the estimated model;}
#'  \item{residuals}{Deviance residuals;}
#'  \item{model.info}{ Model details (equation). The relationship between 
#' the death rate at age x, and the probability of dying between birth and 
#' age 5.}
#' @references \insertAllCited{}
#' @seealso \code{\link{wilmothLT}}
#' @examples 
#' \dontrun{
#' # DATA
#' HMD719f <- HMD719[HMD719$sex == "female", ]
#' 
#' # Fit model
#' x <- c(0,1, seq(5, 110, by = 5))
#' W <- wilmoth(x, LT = HMD719f)
#' }
#' @export  
wilmoth <- function(x, 
                    mx = NULL, 
                    LT = NULL, 
                    verbose = TRUE, 
                    control = list()) {
  
  control  <- do.call("wilmoth_control", control)
  input    <- c(as.list(environment()))
  check_wilmoth_input(I = input)
  info     <- "log m[x] = a[x] + b[x] h + c[x] h^2 + k v[x]"
  na       <- length(x)
  gr_names <- paste0("[", x,",", c(x[-1], "+"), ")")
  
  # Start clock
  if (verbose) {
    N <- if (!is.null(mx)) ncol(mx) else nrow(LT)/na
    pb <- startpb(0, N + 1) 
    on.exit(closepb(pb))
    setpb(pb, 1) 
  }
  
  # Prepare data
  if (!is.null(LT)) {
    mxM <- array(LT$mx, dim = c(na, nrow(LT)/na))
    lxM <- array(LT$lx, dim = c(na, nrow(LT)/na))
    exM <- array(LT$ex, dim = c(na, nrow(LT)/na))
    dimnames(lxM) <- dimnames(mxM) <- dimnames(exM) <- list(x, 1:ncol(mxM))
  }
  
  if (!is.null(mx)) {
    mxM <- as.matrix(mx)
    dimnames(mxM) <- list(x, 1:ncol(mx))
    lxM <- as.matrix(convertFx(x = x, data = mxM, from = "mx", to = "lx"))
    exM <- as.matrix(convertFx(x = x, data = mxM, from = "mx", to = "ex"))
    dimnames(lxM) <- dimnames(exM) <- dimnames(mxM)
  }
  
  # Fit log-quadratic portion of model
  q0_5    <- 1 - lxM["5", ] / lxM["0", ]
  x.f     <- cbind(log(q0_5), log(q0_5)^2)
  y.f     <- t(log(mxM))
  y.f[y.f == -Inf] = -10
  bifit.f <- with(control, bifit(x = x.f, 
                                 y = y.f, 
                                 c = c, 
                                 intercept = intercept, 
                                 tol.biweight = tol.biweight, 
                                 tol.lsfit = tol.lsfit))
  coef    <- data.frame(t(bifit.f$coef))
  
  # Compute residuals and fit SVD portion of model
  y.hat    <- cbind(1, x.f) %*% t(coef)
  resid.f  <- y.hat - y.f
  coef$vx  <- with(control, pmax(0, svd(x = resid.f, nu = nu, nv = nv)$v))
  dimnames(y.hat) <- dimnames(y.f)
  dimnames(coef)  <- list(gr_names, c("ax", "bx", "cx", "vx"))
  
  mx_k <- with(control, 
               find_mx_optim(
                 x = x, 
                 ex = exM, 
                 x.f = x.f, 
                 coef = coef, 
                 radix = radix, 
                 k.int = k.int, 
                 verbose = verbose)
               )
  
  mxfit   <- mx_k$mx
  mxresid <- mxM - mxfit
  
  # Exit
  out <- list(input = input,
              call = match.call(),
              coefficients = coef, 
              k = mx_k$k,
              fitted.values = t(mxfit), 
              residuals = t(mxresid), 
              model.info = info)
  out <- structure(class = 'wilmoth', out)
  return(out)
}


#' Estimate Wilmoth Model Life Table
#' 
#' Construct model life tables based on the Log-Quadratic (wilmoth) estimates
#' with various choices of 2 input parameters: 
#' \code{q0_5, q0_1, q15_45, q15_35, e0} and \code{k}. There are 13 possible 
#' combinations (see examples below).
#' 
#' @details Due to limitations of the R language the notation for probability 
#' of dying \code{nqx} is written \code{qx_n}, where \code{x} and \code{n} are 
#' integers. For example \code{45q15} is represented as \code{q45_15}.
#' 
#' @param object An object of class \code{\link{wilmoth}};
#' @param q0_5 5q0. The probability that a new-born will die during the 
#' subsequent 5 years;
#' @param q0_1 1q0. The probability that a life aged 0 will die during the 
#' following year;
#' @param q15_45 45q15. The probability that a life aged 15 will die during 
#' the subsequent 45 years;
#' @param q15_35 35q15. The probability that a life aged 15 will die during 
#' the subsequent 35 years;
#' @param e0 Life expectancy at birth;
#' @param k k-parameter in the log-quadratic model;
#' @param radix Life table radix. Default: 10^5;
#' @param tol Tolerance level for convergence. The tolerance level, is relevant 
#' for case 12 and 13 (e0 and 45q15 or 35q15 are known);
#' @param maxit Maximum number of iterations allowed. Default: 100;
#' @param ... Additional arguments affecting the predictions produced.
#' @return The output is of class \code{wilmothLT} with the components:
#'  \item{lt}{ Life table matching given inputs}
#'  \item{values}{ Associated values of \code{q0_5, q0_1, q15_45, q15_35, e0} 
#' and \code{k}.}
#' @seealso \code{\link{wilmoth}}
#' @examples 
#' # DATA
#' HMD719f <- HMD719[HMD719$sex == "female", ]
#' hmd_sample <- HMD719f[1:2000, ] 
#' 
#' # Fit Log-quadratic model
#' x <- c(0,1, seq(5, 110, by = 5))
#' W <- wilmoth(x = x, LT = hmd_sample)
#' # use all data to replicate Wilmoth et. al results
#' 
#' # Build life tables with various choices of 2 input parameters
#' 
#' # case 1: Using 5q0 and k
#' L1 <- wilmothLT(W, q0_5 = 0.05, k = 0.1)
#' L1
#' ls(L1)
#'  
#' # case 2: Using 5q0 and e0
#' L2 <- wilmothLT(W, q0_5 = 0.05, e0 = 65)
#' 
#' # case 3: Using 5q0 and 45q15
#' L3 <- wilmothLT(W, q0_5 = 0.05, q15_45 = 0.2)
#' 
#' # case 4: Using 5q0 and 35q15
#' L4 <- wilmothLT(W, q0_5 = 0.05, q15_35 = 0.125)
#' 
#' # case 5: Using 1q0 and k
#' L5 <- wilmothLT(W, q0_1 = 0.01, k = 0.1)
#' 
#' # case 6: Using 1q0 and e0
#' L6 <- wilmothLT(W, q0_1 = 0.01, e0 = 65)
#' 
#' # case 7: Using 1q0 and 45q15
#' L7 <- wilmothLT(W, q0_1 = 0.05, q15_45 = 0.2)
#' 
#' # case 8: Using 1q0 and 35q15
#' L8 <- wilmothLT(W, q0_1 = 0.05, q15_35 = 0.125)
#' 
#' # case 9: Using k and e0
#' L9 <- wilmothLT(W, k = 0.01, e0 = 65)
#' 
#' # case 10: Using k and 45q15
#' L10 <- wilmothLT(W, k = 0.01, q15_45 = 0.2)
#' 
#' # case 11: Using k and 35q15
#' L11 <- wilmothLT(W, k = 0.01, q15_35 = 0.125)
#' 
#' # case 12: Using 45q15 and e0
#' L12 <- wilmothLT(W, q15_45 = 0.125, e0 = 65)
#' 
#' # case 13: Using 35q15 and e0
#' L13 <- wilmothLT(W, q15_35 = 0.15, e0 = 65)
#' 
#' @export
wilmothLT <- function(object, 
                      q0_5 = NULL, 
                      q0_1 = NULL, 
                      q15_45 = NULL, 
                      q15_35 = NULL, 
                      e0 = NULL, 
                      k = NULL, 
                      radix = 1e5, 
                      tol = 1e-9, 
                      maxit = 200, 
                      ...) {
  
  my_case <- find_my_case(q0_5, q0_1, q15_45, q15_35, e0, k)
  cf      <- coef(object)
  x       <- object$input$x
  
  # Cases 1-4:  5q0 is known, plus k, e0, 45q15 or 45q15
  if (my_case == "C1") {
    tmp <- lt_logquad(cf, x, q0_5, k, radix)
  }
  if (my_case %in% c("C2", "C3", "C4")) {
    if (my_case == "C2") fun.k <- function(k) {
      lt_logquad(cf, x, q0_5, k, radix)$lt$ex[1] - e0
    }
    if (my_case == "C3") fun.k <- function(k) { 
      lt <- lt_logquad(cf, x, q0_5, k, radix)$lt
      return(1 - lt[lt$x == 60, "lx"] / lt[lt$x == 15, "lx"] - q15_45)
    }
    if (my_case == "C4") fun.k <- function(k) { 
      lt <- lt_logquad(cf, x, q0_5, k, radix)$lt
      return(1 - lt[lt$x == 50, "lx"] / lt[lt$x == 15, "lx"] - q15_35)
    }
    
    root <- uniroot(f = fun.k, interval = c(-10, 10))$root
    tmp  <- lt_logquad(cf, x, q0_5, k = root, radix) 
  }
  
  # Cases 5-8: 1q0 is known, plus k, e0, 45q15 or 35q15;
  # after finding 5q0 (assume k=0, but it doesn't matter), 
  # these become Cases 1-4
  if (my_case %in% c("C5","C6","C7","C8") ) {
    
    fun.q0_5 <- function(q0_5) {
      lt_logquad(cf, x, q0_5, k = 0, radix)$lt$qx[1] - q0_1
    }
    
    root <- uniroot(f = fun.q0_5, interval = c(1e-5, 0.8))$root
  }
  
  if (my_case == "C5") {
    tmp <- wilmothLT(object, q0_5 = root, k = k, ...)
  }
  
  if (my_case == "C6") {
    tmp <- wilmothLT(object, q0_5 = root, e0 = e0, ...)
  }
  
  if (my_case == "C7") {
    tmp <- wilmothLT(object, q0_5 = root, q15_45 = q15_45, ...)
  }
  
  if (my_case == "C8") {
    tmp <- wilmothLT(object, q0_5 = root, q15_35 = q15_35, ...)
  }
  
  # Cases 9-11: k is known, plus e0, 45q15 or 35q15; 
  # must find 5q0
  if (my_case %in% c("C9", "C10", "C11")) {
    
    if (my_case == "C9") {
      fun.q0_5 = function(q0_5) { 
        lt_logquad(cf, x, q0_5, k, radix)$lt$ex[1] - e0 
      }
    }
    
    if (my_case == "C10") {
      fun.q0_5 = function(q0_5) { 
        lt <- lt_logquad(cf, x, q0_5, k, radix)$lt
        return(1 - lt[lt$x == 60, "lx"] / lt[lt$x == 15, "lx"] - q15_45)
      }
    }
    
    if (my_case == "C11") {
      fun.q0_5 <- function(q0_5) {
        lt <- lt_logquad(cf, x, q0_5, k, radix)$lt
        return(1 - lt[lt$x == 50, "lx"] / lt[lt$x == 15, "lx"] - q15_35)
      }
    } 
    
    root <- uniroot(f = fun.q0_5, interval = c(1e-4, 0.8))$root
    tmp  <- lt_logquad(cf, x, q0_5 = root, k, radix)
  }
  
  # Case 12 and 13: e0 and 45q15 or 35q15 are known; must find both 5q0 and k
  if (my_case %in% c("C12", "C13")) {
    k    <- q0_5 <- 0
    iter <- crit <- 1
    
    while (crit > tol & iter <= maxit) {
      k.old    <- k
      q0_5.old <- q0_5
      # Get new 5q0 from e0 assuming k
      q0_5 <- wilmothLT(object, k = k, e0 = e0, ...)$values$q0_5 
      # Get k from 45q15 or 35q15 assuming 5q0
      if (my_case == "C12") {
        tmp = wilmothLT(object, q0_5 = q0_5, q15_45 = q15_45, ...)
      }
      
      if (my_case == "C13") {
        tmp = wilmothLT(object, q0_5 = q0_5, q15_35 = q15_35, ...)
      }
      
      k    <- tmp$values$k
      crit <- sum(abs(c(k, q0_5) - c(k.old, q0_5.old)))
      iter <- iter + 1
    }
    
    if (iter > maxit) {
      warning("number of iterations reached maximum without convergence", 
              call. = FALSE)
    }
  }
  
  # Return life table plus values of the 6 possible inputs
  out = list(lt = tmp$lt, 
             values = tmp$values)
  out = structure(class = "wilmothLT", out)
  return(out)
}



#' Compute WLS using Tukey's biweight function
#' 
#' @param x A matrix whose rows correspond to cases and whose columns 
#' correspond to variables. x matrix should not include a column of ones
#' @param y The responses, possibly a matrix if you want to fit multiple left 
#' hand sides. y can be matrix if there are multiple left-hand sides
#' @param c Tuning constant when compute weighted least squared using 
#' Tukey's biweight function. Typically, between 6 to 9.
#' @param intercept logical. Whether or not an intercept term should be used in 
#' least squares estimation. See \code{\link{lsfit}}.
#' @param tol.biweight Tolerance for convergence of the biweight fit
#' @param tol.lsfit Tolerance to be used in the matrix decomposition in 
#' least squares estimation. See \code{\link{lsfit}}.
#' @keywords internal
bifit <- function(x, 
                  y, 
                  c, 
                  intercept, 
                  tol.biweight, 
                  tol.lsfit) {
  
  if (is.vector(y)) {
    coef.old <- 0
    iter <- 0
    coef <- 1 
    wt   <- NULL
    
    while (max(abs(coef - coef.old)) > tol.biweight) {
      iter <- iter + 1
      z <- lsfit(x = x, 
                 y = y, 
                 wt = wt, 
                 intercept = intercept, 
                 tolerance = tol.lsfit)
      
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
      z <- bifit(x = x, 
                 y = y[,j], 
                 c = c, 
                 intercept = intercept, 
                 tol.biweight = tol.biweight, 
                 tol.lsfit = tol.lsfit)
      
      resid <- cbind(resid, z$resid)
      wt    <- cbind(wt, z$wt)
      u     <- cbind(u, z$u)
      coef  <- cbind(coef, z$coef) 
    }
    
    dimnames(resid) <- dimnames(wt) <- dimnames(u) <- dimnames(y)
    dimnames(coef)  <- list(
      paste0("b", 0:(ncol(cbind(1, x)) - 1)), dimnames(y)[[2]]
      )
    
  }
  out <- list(coef = coef, 
              residuals = resid, 
              intercept = intercept, 
              wt = wt, 
              u = u)
  return(out)
}


#' Compute fitted mx
#' 
#' @inheritParams wilmoth
#' @inheritParams wilmothLT
#' @param ex Life expectancy matrix
#' @param x.f Data frame containing the q0_5 values.
#' @param coef Estimated coefficients (ax, bx, cx) of the Log-Quadratic model.
#' @param k.int a vector containing the end-points of the interval to be 
#' searched for k parameter.
#' @keywords internal
find_mx_optim <- function(x, 
                          ex, 
                          x.f, 
                          coef, 
                          radix, 
                          k.int, 
                          verbose) {
  
  fit <- ex * 0
  k   <- rep(NA, nrow(x.f))
  N   <- length(k)
  if (verbose) pb <- startpb(1, N + 1)
  
  for (i in 1:N) {
    if (verbose) setpb(pb, i + 1)
    q0_5     <- exp(x.f[i, 1])
    e0t      <- ex[1, i]
    fun.k    <- function(k) {
      round(lt_logquad(coef, x, q0_5, k, radix)$values$e0 - e0t, 4)
    }
    k[i]     <- uniroot(f = fun.k, interval = k.int)$root
    fit[, i] <- lt_logquad(coef, x, q0_5, k = k[i], radix)$lt$mx
  }
  
  # Exit
  out <- list(mx = fit, k = k)
  return(out)
}


#' Estimated life table using the log-quadratic model
#' 
#' @param coefs Estimated coefficients
#' @inheritParams wilmoth
#' @inheritParams wilmothLT
#' @keywords internal
#' @export
lt_logquad <- function(coefs, 
                          x, 
                          q0_5, 
                          k, 
                          radix) {
  
  h     <- log(q0_5)
  mx    <- with(as.list(coefs), exp(ax + bx*h + cx*h^2 + vx*k))
  # Force 4q1 (and thus 4m1) to be consistent with 1q0 and 5q0
  qx    <- mx_qx(x, mx, out = "qx")
  qx[2] <- 1 - (1 - q0_5)/(1 - qx[1])
  mx[2] <- mx_qx(x, qx, out = "mx")[2]
  names(mx) = names(qx) <- rownames(coefs)
  
  LT     <- LifeTable(x = x, mx = mx, lx0 = radix)$lt
  e0     <- LT$ex[1]
  q0_1   <- LT$qx[1]
  q15_45 <- 1 - LT[LT$x == 60, "lx"] / LT[LT$x == 15, "lx"]
  q15_35 <- 1 - LT[LT$x == 50, "lx"] / LT[LT$x == 15, "lx"]
  values <- data.frame(k, q0_1, q0_5, q15_35, q15_45, e0, row.names = "")
  
  # Exit
  out <- list(lt = LT, values = values)
  return(out)
}


#' Check wilmoth arguments
#' 
#' @param input a list containing the input objects of the wilmoth function
#' @keywords internal
check_wilmoth_input <- function(I) {
  
  if (!is.numeric(I$x)) {
    stop("'x' must be of class 'numeric'", call. = FALSE)
  }
  
  if (!is.null(I$mx)) {
    if (length(I$x) != nrow(I$mx)) {
      stop("length(x) must be identical with nrow(mx)", call. = FALSE)
    }
  }  
}

#' Function that determines the case/problem we have to solve
#' It also performs some checks
#' @inheritParams wilmothLT
#' @keywords internal
find_my_case <- function(q0_5, q0_1, q15_45, q15_35, e0, k) {
  # Test that at least of one of 1q0 and 5q0 is null
  input   <- as.list(environment())
  my_case <- unlist(lapply(input, is.null))
  
  if (sum(my_case[c(1, 2)]) == 0) {
    stop("cannot have both 'q0_1' and 'q0_5' as inputs", call. = FALSE)
  }
  
  # Test that at least of one of 45q15 and 35q15 is null
  if (sum(my_case[c(3, 4)]) == 0) {
    stop("cannot have both 'q15_45' and 'q15_35' as inputs", call. = FALSE)
  }
  
  # Test that exactly two inputs are non-null
  if (sum(my_case) != 4) {
    stop("must have exactly two inputs", call. = FALSE)
  }
  
  # There are 13 cases:  "6 choose 2" = 15, but we disallow two cases 
  # (1q0 and 5q0, or 45q15 and 35q15)
  if (all(my_case == c(F,T,T,T,T,F))) case = "C1" 
  if (all(my_case == c(F,T,T,T,F,T))) case = "C2"
  if (all(my_case == c(F,T,F,T,T,T))) case = "C3" 
  if (all(my_case == c(F,T,T,F,T,T))) case = "C4"
  if (all(my_case == c(T,F,T,T,T,F))) case = "C5" 
  if (all(my_case == c(T,F,T,T,F,T))) case = "C6"   
  if (all(my_case == c(T,F,F,T,T,T))) case = "C7"
  if (all(my_case == c(T,F,T,F,T,T))) case = "C8" 
  if (all(my_case == c(T,T,T,T,F,F))) case = "C9" 
  if (all(my_case == c(T,T,F,T,T,F))) case = "C10" 
  if (all(my_case == c(T,T,T,F,T,F))) case = "C11" 
  if (all(my_case == c(T,T,F,T,F,T))) case = "C12" 
  if (all(my_case == c(T,T,T,F,F,T))) case = "C13" 
  return(case)
}

#' Auxiliary for Controlling \code{wilmoth} Fitting.
#' 
#' @inheritParams bifit
#' @inheritParams find_mx_optim
#' @param nu The number of left singular vectors to be computed in SVD method. 
#' This must between 0 and n = nrow(x). See \code{\link{svd}}.
#' @param nv the number of right singular vectors to be computed. 
#' This must be between 0 and p = ncol(x). See \code{\link{svd}}.
#' @return List with control parameters.
#' @seealso \code{\link{wilmoth}}
#' @export
wilmoth_control <- function(c = 6, 
                            intercept = TRUE, 
                            tol.biweight = 1e-06, 
                            tol.lsfit = 1e-07,
                            radix = 1e5, 
                            k.int = c(-25, 20),
                            nu = 1,
                            nv = 1) {
  
  out <- c(as.list(environment()))
  if (tol.biweight <= 0) {
    stop("'tol.biweight' should be greater than 0")
  }
  
  return(out)
}



# ----------------------------------------------
# S3 functions

#' Print function for wilmoth
#' 
#' @param x An object of class \code{wilmoth}.
#' @param ... Further arguments passed to or from other methods. 
#' @keywords internal
#' @export
print.wilmoth <- function(x, ...){
  cat("\nThe Log-Quadratic model:\n")
  cat(x$model.info, "\n")
  cat("\nCoefficients:\n")
  print(round(coef(x), 4))
}

#' Summary for wilmoth
#' 
#' @param object An object of class \code{wilmoth}.
#' @inheritParams print.wilmoth
#' @keywords internal
#' @export
summary.wilmoth <- function(object, ...) {
  print.wilmoth(x = object, ...)
}



