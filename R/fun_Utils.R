

#' Summary function - display head and tail in a single data.frame
#' The code for this function was first written for 'psych' R package
#' @importFrom utils head tail
#' @keywords internal
head_tail <- function(x, hlength = 4, tlength = 4, digits = 2, ellipsis = TRUE) 
{
  if (is.data.frame(x) | is.matrix(x)) {
    if (is.matrix(x)) 
      x <- data.frame(unclass(x))
    nvar <- dim(x)[2]
    dots <- rep("...", nvar)
    h <- data.frame(head(x, hlength))
    t <- data.frame(tail(x, tlength))
    for (i in 1:nvar) {
      if (is.numeric(h[1, i])) {
        h[i] <- round(h[i], digits)
        t[i] <- round(t[i], digits)
      }
      else {
        dots[i] <- NA
      }
    }
    if (ellipsis) {
      head.tail <- rbind(h, ... = dots, t)
    }
    else {
      head.tail <- rbind(h, t)
    }
  }
  else {
    h <- head(x, hlength)
    t <- tail(x, tlength)
    if (ellipsis) {
      head.tail <- rbind(h, "...       ...", t)
    }
    else {
      head.tail <- rbind(h, t)
      head.tail <- as.matrix(head.tail)
    }
  }
  return(head.tail)
}

