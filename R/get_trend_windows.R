#' Consensus ADMM for overlapping windows
#'
#' \code{consensus ADMM} Uses ADMM to smooth overlapping windows
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param lambda smoothing penalty parameter
#' @param k order of differencing
#' @param window_size size of windows to use
#' @param overlap integer length of overlap between windows
#' @param max_iter Maximum number of iterations
#' @param rho parameter for ADMM
#' @param update number of iterations at which to print residuals
#' @param use_gurobi TRUE if gurobi solver is installed and should be used
#' @param eps_abs absolute threshold for stopping criteria
#' @param eps_rel relative threshold for stopping criteria
#' @examples
#' require(Matrix)
#' n <- 100
#' x <- seq(1, n, 1)
#' y <- sin(x*2*pi/n) + rnorm(n, 0, .4)
#' lambda <- 10
#' k <- 3
#' y_n <- length(y)
#' overlap <- 20
#' window_size <- round((y_n+overlap)/2)
#' tau <- c(0.05, .2)
#' max_iter <- 20
#' trend <- get_trend_windows(y, tau, lambda, k, window_size, overlap, max_iter)
#' plot(trend[,1]~x, type="l")
#' @export
get_trend_windows <- function(y, tau, lambda, k, window_size,
                           overlap, max_iter, rho=1, update=10, 
                           use_gurobi = TRUE, 
                           eps_abs = .05, 
                           eps_rel = 1e-3){
  min_y <- min(y, na.rm=T)
  max_y <- max(y, na.rm=T)
  y <- 200*(y-min_y)/(max_y-min_y)
  if (use_gurobi){
    solver <- "gurobi"
  } else {
    # First estimate uses LP not QP
    solver <- "lpSolve"
  }
  
  window_size <- round(window_size)
  
  y_n <- length(y)
  tau <- sort(tau)
  nT <- length(tau)
  D <- get_Dk(window_size, k)
  n_windows <- ceiling((y_n-overlap)/(window_size-overlap))
  windows <- matrix(FALSE, y_n, n_windows)
  y_list <- list()
  w_list <- list()
  
  # Initial values
  for (i in 1:n_windows){
    start_I <- 1+(window_size-overlap)*(i-1)
    end_I <- min((window_size + (window_size-overlap)*(i-1)), y_n)
    windows[start_I:end_I,i] <- TRUE
    y_list[[i]] <- y[start_I:end_I]
    len <- end_I - start_I + 1
    w_list[[i]] <- rep(0, (2*len - k)*length(tau))
  }
  overlapInd <- rowSums(windows) > 1
  
  # Window initial LP fit
  model_list <- lapply(y_list, get_model, tau=tau, lambda=lambda, k=k)
  phi_list <- mapply(solve_model, model_list, y_list, 
                     solver = solver, trend=FALSE, SIMPLIFY = FALSE)
  eta_list <- mapply(get_eta, phi_list, y_list, k=k, SIMPLIFY = FALSE)
  
  # Change to QP solver
  if (solver == "lpSolve"){
    solver <- "quadprog"
  }
  # Consensus update
  phiBar_list <- update_consensus(phi_list, windows, overlapInd)
  etaBar_list <- mapply(get_eta, phiBar_list, y_list, k=k, SIMPLIFY = FALSE)
  
  # Use w_list for z since it is all zeros, don't want to store current z
  eta0 <- etaBar_list
  phi0 <- phiBar_list
  for (i in 1:length(eta0)){
    eta0[[i]][] <- 0
    phi0[[i]][] <- 0
  }
  model_list <- mapply(update_model, model_list, w_list, phi0,
                       eta0, rho=rho, nT=nT, SIMPLIFY = FALSE)

  # Dual update
  w_list <- mapply(update_dual, w_list, 
                   phi_list, phiBar_list, 
                   eta_list, etaBar_list, 
                   MoreArgs = list(rho=rho), SIMPLIFY = FALSE)

  dual_norm <- double(max_iter)
  primal_norm <- double(max_iter)
  phiBar_listk <- phiBar_list
  iter <- 1

  while(iter <= max_iter){

    # Window update
    phi_list <- update_windows(w_list, phiBar_list, etaBar_list,
                               model_list, rho, nT, solver)
    
    # Consensus update
    phiBar_list <- update_consensus(phi_list, windows, overlapInd)
    etaBar_list <- mapply(get_eta, phiBar_list, y_list, k=k, SIMPLIFY = FALSE)
    
    # Dual update
    w_list <- mapply(update_dual, w_list, 
                     phi_list, phiBar_list, 
                     eta_list, etaBar_list, 
                     MoreArgs = list(rho=rho), SIMPLIFY = FALSE)
    
    # Convergence Metrics

    dual_norm[iter] <- rho*sqrt(sum(mapply(list_diff_norm, phiBar_list, phiBar_listk)))
    primal_norm[iter] <- sqrt(sum(mapply(list_diff_norm, phi_list, phiBar_list)))
    eps_pri <- sqrt(nT*y_n)*eps_abs + 
      eps_rel*max(c(sapply(phi_list, Matrix::norm, type="F"), 
                    sapply(phiBar_list, Matrix::norm, type="F")))
    
    eps_dual <- sqrt(nT*y_n)*eps_abs + 
      eps_rel*sqrt(sum(sapply(w_list, norm, type="2")^2))
    
    phiBar_listk <- phiBar_list
    
    if (iter %% update == 0){
      print(sprintf("Iteration: %d Primal Resid Norm: %.4f eps_pri: %.4f, Dual Resid Norm: %.4f  eps_dual %.4f", 
                    iter,
                    primal_norm[iter],
                    eps_pri,
                    dual_norm[iter],
                    eps_dual))

    }

    if (iter == 1 ){
      primal_resid0 <-  primal_norm[iter] 
      dual_resid0 <- dual_norm[iter]
    } else if(dual_norm[iter] < eps_dual & 
              primal_norm[iter] < eps_pri){
      primal_norm <- primal_norm[1:iter]
      dual_norm <- dual_norm[1:iter]
      print(sprintf("Converged in %d iterations", iter))
      break
    }
    iter <- iter+1
  }
  y[is.na(y)] <- 0
  theta <- ((y - get_phiBar(phiBar_list, windows))/200)*(max_y-min_y) + min_y
  return(theta)
}

list_diff_norm <- function(list1, list2) {
  return(Matrix::norm(list1 - list2, type = "F"))
}
