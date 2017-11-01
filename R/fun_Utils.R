# Function borrowed from MortalityLaws R package that are not exported yet.


#' Convert 'mx' into 'qx' and viceversa.
#'
#' Function to convert mx into qx and back, using the constant force of 
#' mortality assumption (CFM).
#' @param ux a vector of mx or qx
#' @keywords internal
mx_qx <- function(x, ux, out = "qx"){
  if (!(out %in% c("qx", "mx"))) stop("out must be: 'qx' or 'mx'", call. = FALSE)
  N     <- length(x)
  nx    <- c(diff(x), Inf)
  if (out == "qx") {
    eta = 1 - exp(-nx*ux)
    eta[is.na(ux)] <- 1
    eta[x >= 100 & ux == 0]  <- 1
    if (max(x) > 100) eta[N] <- 1
  }
  if (out == "mx") {
    eta = -log(1 - ux)/nx
    eta[is.infinite(eta)] <- max(eta[!is.infinite(eta)], na.rm = T)
    eta[is.na(eta)] <- max(eta, na.rm = T)
    # here if qx[N] = 1 then mx[N] = NaN therefore we apply a simple extrapolation method
    eta[N] = eta[N - 1]^2 / eta[N - 2]
  }
  return(eta)
}


#' Summary function - display head and tail in a single data.frame
#' @param x A matrix or data frame or free text
#' @param hlength The number of lines at the beginning to show
#' @param tlength The number of lines at the end to show
#' @param digits Round off the data to digits
#' @param ellipsis separate the head and tail with dots
#' @keywords internal
head_tail <- function(x, hlength = 4, tlength = 4, digits = 4, ellipsis = TRUE){
  if (is.data.frame(x) | is.matrix(x)) {
    if (is.matrix(x)) x = data.frame(unclass(x))
    nvar <- dim(x)[2]
    dots <- rep("...", nvar)
    h    <- data.frame(head(x, hlength))
    t    <- data.frame(tail(x, tlength))
    for (i in 1:nvar) {
      if (is.numeric(h[1, i])) {
        h[i] <- round(h[i], digits)
        t[i] <- round(t[i], digits)
      } else {
        dots[i] <- NA
      }
    }
    out <- if (ellipsis) rbind(h, ... = dots, t) else rbind(h, t)
  } else {
    h <- head(x, hlength)
    t <- tail(x, tlength)
    out <- paste(paste(h, collapse = " "), "...   ...", paste(t, collapse = " "))
  }
  return(out)
}

