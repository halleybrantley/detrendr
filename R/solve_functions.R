#' Solves the model using chosen solver
#'
#' \code{solve_model} Returns the trend
#'
#' @param model  model object containing coefficients and constraints
#' @param solver solver to be used, options are "gurobi", "Rglpk", and "lpSolve"
#' @param y vector of observations, only needed if trend = TRUE
#' @param trend if TRUE returns trend, if FALSE returns residuals
#' @export
solve_model <- function(model, solver, y=NULL, trend = TRUE){
  if (trend & is.null(y)){
    stop("y required to get trend")
  }
  
  if (solver == "gurobi"){
    x <- solve_gurobi(model)
  } else if (solver == "Rglpk") {
    x <- solve_glpk(model)
  } else if (solver == "lpSolve"){
    x <- solve_lp(model)
  } else if (solver == "ipop"){
    x <- solve_ipop(model)
  } else if (solver == "quadprog"){
    x <- solve_quadprog(model)
  } else {
    stop("Solver must be one of 'gurobi', 'Rglpk', 'lpSolve', 'ipop', 'quadprog'")
  }
  
  nT <- model$nT
  n <- model$n
  np <- model$np
  phi <- matrix(0, nrow=n, ncol = nT)
  
  if (!is.null(x)){
    for (i in 1:nT){
      phi[,i] <- x[(1+np*(i-1)):(n+np*(i-1))] - 
        x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
    }
  }

  if (trend){
    y[which(is.na(y))] <- 0
    return(y-phi)
  } else {
    return(phi)
  }
}  


# Wrappers for different solvers

solve_gurobi <- function(model){
  params <- list(OutputFlag=0)
  result <- gurobi::gurobi(model, params)
  if (result$status == "NUMERIC"){
    print("Optimization was terminated due to unrecoverable numerical difficulties.")
    result <- gurobi::gurobi(model, params)
  }
  return(result$x)
}

solve_glpk <- function(model_list){
  result <- Rglpk::Rglpk_solve_LP(obj = model_list$obj, 
                           mat = model_list$A, 
                           dir = paste0(model_list$sense, "="), 
                           rhs = model_list$rhs)
  return(result$solution)
}

solve_lp <- function(model_list){
  result <- lpSolve::lp(objective.in = model_list$obj, 
               const.mat = as.matrix(model_list$A), 
               const.dir = model_list$sense, 
               const.rhs = model_list$rhs)
  return(result$solution)
}

solve_ipop <- function(model){
  
  b <- model$rhs
  r <- as.numeric(sub(">", 100, sub("=", 1e-8, model$sense)))
  H <- as.matrix(2*model$Q)
  result <- kernlab::ipop(
    H = H,
    c = model$obj,
    A = as.matrix(model$A), 
    b = b, 
    r = r,
    l = rep(0, length(model$obj)),
    u = rep(1000, length(model$obj)), 
    sigf = 1e-20, 
    margin = 0.05
  )
  output <- kernlab::primal(result)
  return(output)
}

solve_quadprog <- function(model){
  Amat <- Matrix::t(rbind(model$A, Matrix::Diagonal(length(model$obj))))
  bvec <- c(model$rhs, rep(0, length(model$obj)))
  result <- quadprog::solve.QP(
    Dmat = 2*model$Q, 
    dvec = -model$obj,
    Amat = Amat, 
    bvec = bvec,
    meq = sum(model$sense == "=")
  )
  return(result$solution)
}

