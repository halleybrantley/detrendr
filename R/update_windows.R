#' Update windows
#'
#' \code{update_windows} Minimizes the Lagrangian for each window
#'
#' @param w_list Dual variable
#' @param phiBar_list Consensus variable for phi
#' @param etaBar_list Consensus variable for eta
#' @param model_list List of model objects for solver
#' @param rho ADMM parameter
#' @param nT number of quantiles
#' @param solver must be "gurobi" or "quadprog"
#' @export
update_windows <- function(w_list, phiBar_list, etaBar_list, 
                           model_list, rho, nT, solver="gurobi"){
 
  model_list <- mapply(update_model, model_list, w_list, phiBar_list,
                       etaBar_list,
                       rho=rho, nT=nT, SIMPLIFY = FALSE)

  phi_list <- lapply(model_list, solve_model, solver=solver, trend = FALSE)
  # If solver didn't succeed use value of phiBar
  phi_list0 <- mapply(function(x,y) {if(!any(x != 0)){return(y)}else{return(x)}}, 
                     phi_list, phiBar_list)
  return(phi_list)
}

#' Update consensus variable (overlap)
#'
#' \code{update_consensus}
#'
#' @param phi_list Primal variables
#' @param windows Matrix indicating window group
#' @param overlapInd Indices of overlap between windows
#' @export
update_consensus <- function(phi_list, windows, overlapInd){
  nT <- dim(phi_list[[1]])[2]
  phiBar <- matrix(0, nrow = nrow(windows), ncol = nT)
  n_windows <- length(phi_list)
  
  for (i in 1:n_windows){
    phiBar[windows[,i],] <- phiBar[windows[,i],] + phi_list[[i]]
  }
  phiBar[overlapInd, ] <- phiBar[overlapInd, ]/2

  phiBar_list <- list()
  for (i in 1:n_windows){
    phiBar_list[[i]] <- phiBar[windows[,i],,drop=FALSE]
  }
  return(phiBar_list)
}

#' Update consensus variable (overlap)
#'
#' \code{get_phiBar}
#'
#' @param phiBar_list List of consensus variable
#' @param windows Matrix indicating window group
#' @export
get_phiBar <- function(phiBar_list, windows){
  phiBar <- matrix(0, nrow = nrow(windows), ncol = ncol(phiBar_list[[1]]))
  for (i in 1:length(phiBar_list)){
     phiBar[windows[,i],] <- phiBar_list[[i]] 
  }
  return(phiBar)
}


#' Update dual variable (overlap)
#'
#' \code{update_dual}
#'
#' @param w dual variable
#' @param phi primal variable
#' @param phiBar consensus variable
#' @param eta primal variable equal to Dtheta
#' @param etaBar consensus variable for eta
#' @param rho ADMM step size parameter
#' @export
update_dual <- function(w, phi, phiBar, eta, etaBar, rho) {
  w + as.numeric(rbind(rho*(phi - phiBar), rho*(eta-etaBar)))
}

#' Update model with ADMM step parameter and dual and consensus variables
#'
#' \code{update_model}
#'
#' @param model List of QP model elements
#' @param w dual variable
#' @param phiBar consensus variable for phi
#' @param etaBar consensus variable for eta
#' @param rho ADMM step size
#' @param nT Number of quantiles being estimated
#' @export
update_model <- function(model, w, phiBar, etaBar, rho, nT){
  z <- as.numeric(rbind(phiBar, etaBar))
  n <- length(phiBar)/nT
  np <- length(model$obj)/nT
  
  for (i in 1:nT){
    objInd <- (1+np*(i-1)):(np*i)
    dual_theta <- (1 + np/2*(i-1)):(np/2*(i-1) + n)
    dual_eta <- (1 + np/2*(i-1) + n):(np/2*i)
    # Add in dual/consensus variables
    model$obj[objInd] <- model$obj[objInd]  +
      c(w[dual_theta]-rho*z[dual_theta], -(w[dual_theta]-rho*z[dual_theta]),
        w[dual_eta]-rho*z[dual_eta], -(w[dual_eta]-rho*z[dual_eta]))
  }
  
  if (is.null(model$Q)) {
    model$Q <- Matrix::Diagonal(length(model$obj), rho/2)
  }
  return(model)
}

#' Update dual variable (overlap)
#'
#' \code{get_eta}
#'
#' @param phi matrix of y-theta
#' @param y vector of observations
#' @param k order of differencing matrix
#' @export
get_eta <- function(phi, y, k){
  D <- get_Dk(length(y), k)
  eta <- D%*%y - D%*%phi
  eta[is.na(eta)] <- 0
  return(as.matrix(eta))
}


