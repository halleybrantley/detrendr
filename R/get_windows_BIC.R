#' Choose smoothing parameter using SIC
#'
#' \code{get_windows_BIC} Selects smoothing parameter using Schwarz Information
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
#' @param window_size integer of size of windows, last window may be shorter
#' @param overlap number of observations in overlap between windows
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @param df_tol tolerance for determining degrees of freedom (Dtheta > df_tol)
#' @param gamma parameter for eBIC
#' @param plot_lambda TRUE/FALSE for plotting lambda by model criteria
#' @param solver LP solver, can be "gurobi" or "lpSolve"
#' @param criteria criteria to use for lambda selection, must be "eBIC", "SIC", or "valid"  
#' @param max_iter Maximum number of iteration
#' @param eps_abs absolute threshold for stopping criteria
#' @param eps_rel relative threshold for stopping criteria
#' @param update number of iterations to print residual values
#' @param rho ADMM parameter
#' @importFrom graphics abline lines plot
#' @importFrom utils install.packages
#' @export
get_windows_BIC <- function(y, tau, k, window_size, overlap,
                             lambdaSeq = exp(seq(0, 14, 1)), 
                             df_tol = 1e-9, 
                             gamma = 1,
                             plot_lambda = FALSE, 
                             solver = NULL, 
                             criteria = "eBIC", 
                             max_iter = 10, 
                             eps_abs = 0.05, 
                             eps_rel = 1e-3, 
                             update = 10, 
                             rho = 5){
  
  if(!(criteria %in% c("eBIC", "valid", "SIC"))){
    stop("criteria must be one of 'eBIC', 'valid', 'SIC'")
  }
  
  # Set linear program solver
  if (is.null(solver)){
    pkgs <- installed.packages()[,"Package"]
    if("gurobi" %in% pkgs){
      use_gurobi <- TRUE
    } else {
      use_gurobi <- FALSE
    }
  } else if (solver == "gurobi"){
    use_gurobi <- TRUE
  } else {
    use_gurobi <- FALSE
  }
  
  # min_y <- min(y, na.rm=T)
  # max_y <- max(y, na.rm=T)
  # y <- 10*(y-min_y)/(max_y-min_y)
  n <- length(y)
  
  n <- length(y)
  m <- n-k
  missInd <- which(is.na(y))  
  
  if (criteria == "valid"){
    validID <- seq(5, length(y), 5)
    yValid <- y[validID]
    y[validID] <- NA
    D <- NULL
  } else {
    validID <- NULL
    yValid <- NULL
    D <- get_Dk(n, k)
  }
  
  df <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
  BIC <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))

  for (i in 1:length(lambdaSeq)){
    f_trend <- get_trend_windows(y, tau, lambdaSeq[i], k=k, window_size = window_size,
                                 overlap=overlap, max_iter=max_iter, update=update, 
                                 use_gurobi = use_gurobi, 
                                 eps_abs = eps_abs, eps_rel = eps_rel, 
                                 rho = rho)
    model_crit <- get_criteria(criteria, f_trend, y, tau, 
                               D, df_tol, gamma, 
                               validID, yValid)
    df[i,] <- model_crit$df
    BIC[i,] <- model_crit$BIC
    print(sprintf("i=%d lambda=%f", i, lambdaSeq[i]))
  }
  
  lambda <- lambdaSeq[apply(BIC, 2, which.min)]
  
  if (plot_lambda){
    plot(BIC[,1]~log(lambdaSeq), type="l", col="red", 
         ylim = c(min(BIC), max(BIC)))
    for (i in 1:length(tau)){
      lines(BIC[,i]~log(lambdaSeq))
    }
    abline(v=log(lambda))
  }
  
  f_trend <- get_trend_windows(y, tau, lambda, k, window_size,
                               overlap, max_iter=max_iter, update=update,
                               eps_abs = eps_abs, eps_rel = eps_rel,
                               rho = rho)
  
  return(list(trend = f_trend, #*(max_y-min_y)/10 + min_y,
              lambda = lambda, 
              BIC = BIC, 
              df = df))
}
