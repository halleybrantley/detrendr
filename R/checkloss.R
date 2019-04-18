# Functions for choosing smoothing parameter

#' Evaluate checkloss function
#'
#' \code{checkloss}
#'
#' @param e argument of checkloss function
#' @param tau quantile to be used
#' @export
checkloss <- function(e, tau){
  if (ncol(e) != length(tau)){
    stop("Number of columns in y must be same as length of tau")
  }
  obj <- e
  for (i in 1:length(tau)){
    obj[,i] <- obj[,i]*tau[i]
    obj[e[,i] < 0,i] <- e[e[,i] < 0,i]*(tau[i]-1)
  }
  return(obj)
  
}