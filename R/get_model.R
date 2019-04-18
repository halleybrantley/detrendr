#' Get objective linear program 
#'
#' \code{get_objective} Returns objective coefficent vector
#'
#' @param tau quantiles being estimated
#' @param lambda smoothing parameters, should be same length as tau
#' @param n length of data
#' @param m n-k where k is order of differencing matrix
#' @param missInd output of which(is.na(y))
#' @export
get_objective <- function(tau, lambda, n, m, missInd){
  obj <- c()
  nT <- length(tau)
  for (i in 1:nT){
    obj0 <- c(rep(tau[i], n), rep((1-tau[i]), n), rep(lambda[i], 2*m))
    obj0[missInd] <- 0
    obj0[missInd + n] <- 0
    obj <- c(obj, obj0)
  }

  return(obj)
}

#' Get constraint matrix for linear program 
#'
#' \code{get_constraint_mat} Returns the constraint matrix
#'
#' @param D discrete difference matrix
#' @param nT number of quantiles being estimated
#' @export
get_constraint_mat <- function(D, nT){
  n <- ncol(D)
  m <- nrow(D)
  np <- 2*n + 2*m
  # Constraint Matrix
  if (nT == 1){
    A  <- cbind(D, -D, Matrix::Diagonal(m), -Matrix::Diagonal(m))
  } else {
    A <- Matrix::Matrix(0, nrow =  m*nT + n*(nT-1), ncol= np*nT, sparse=TRUE)
    for (i in 1:nT){
      # D%*%theta = eta constraint
      A[(1+m*(i-1)):(m*i), (1+np*(i-1)):(np*i)] <-
        cbind(D, -D, Matrix::Diagonal(m), -Matrix::Diagonal(m))
      
      # Non-crossing theta(tau) constrains
      if (i < nT){
        A[(m*nT+1+n*(i-1)):(m*nT+n*(i)), (1+(i-1)*np):(2*n + (i-1)*np)] <-
          cbind(Matrix::Diagonal(n), -Matrix::Diagonal(n))
        
        A[(m*nT+1+n*(i-1)):(m*nT+n*(i)),
          (1+i*np):(2*n + i*np)] <-
          cbind(-Matrix::Diagonal(n), Matrix::Diagonal(n))
      }
    }
  }
  return(A)
}

#' Get linear program arguments for gurobi or glpk solver
#'
#' \code{get_model} Returns the linear program solver arguments
#'
#' @param y observed data, should be equally spaced, may contain NA
#' @param tau vector of quantile levels at which to evaluate trend
#' @param lambda vector of penalty parameters controlling smoothness
#' @param k order of differencing
#' @export
get_model <- function(y, tau, lambda, k){
  
  if(tau >= 1 || tau <= 0){
    stop("tau must be between 0 and 1.")
  }
  if (length(lambda) == 1 && length(tau) != 1) {
    lambda <- rep(lambda, length(tau))
    message("Using same lambda for all quantiles")
  } else if (length(lambda) != length(tau)){
    stop("lambda must be the same size as tau")
  }
  
  tau <- sort(tau)
  D <- get_Dk(length(y), k)
  n <- length(y)
  m <- nrow(D)
  nT <- length(tau)
  missInd <- which(is.na(y))
  y[missInd] <- 0
  
  obj <- get_objective(tau, lambda, n, m, missInd)
  rhs <- c(rep(as.numeric(D%*%y), nT), rep(0, n*(nT-1)))
  A <- get_constraint_mat(D, nT)
  
  # Constraint Types
  sense <- c(rep('=', m*nT), rep(">", n*(nT-1)))
  modelsense <- "min"
  return(list(obj=obj, A=A, rhs=rhs, sense=sense, modelsense = modelsense, 
              n=n, np=2*n+2*m, nT=nT))
}