#' Choose smoothing parameter using SIC
#'
#' \code{get_trend_cv} Selects smoothing parameter using Schwarz Information
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @param df_tol tolerance for determining degrees of freedom (Dtheta > df_tol)
#' @param gamma parameter for eBIC
#' @param plot_lambda TRUE/FALSE for plotting lambda by model criteria
#' @param solver LP solver, can be "gurobi", "Rglpk", or "lpSolve"
#' @param criteria criteria to use for lambda selection, must be "eBIC", "SIC", or 
#' "valid"
#' @examples 
#' require(Matrix)
#' n <- 100
#' x <- seq(1, n, 1)
#' y <- sin(x*2*pi/n) + rnorm(n, 0, .4)
#' lambdaSeq <- exp(seq(-2, 5, 1))
#' k <- 3
#' tau <- c(0.05, .2)
#' trend <- get_trend_BIC(y, tau, k, lambdaSeq, plot_lambda = TRUE)
#' plot(trend$theta[,1]~x, type="l")  
#' @export
get_trend_BIC <- function(y, tau, k,
                       lambdaSeq = exp(seq(0, 14, 1)), 
                       df_tol = 1e-9, 
                       gamma = 1,
                       plot_lambda = FALSE, 
                       solver = NULL, 
                       criteria = "eBIC"){
  n <- length(y)
  if(!(criteria %in% c("eBIC", "valid", "SIC"))){
    stop("criteria must be one of 'eBIC', 'valid', 'SIC'")
  }
  
  # Set linear program solver
  if (is.null(solver)){
    pkgs <- installed.packages()[,"Package"]
    if("gurobi" %in% pkgs){
      solver <- "gurobi"
    } else if ("Rglpk" %in% pkgs){
      solver <- "Rglpk"
    } else {
      solver <- "lpSove"
    }
  }
  
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

  model <- get_model(y, tau, lambdaSeq[1], k)
  m <- n-k
  missInd <- which(is.na(y))
  
  for (i in 1:length(lambdaSeq)){
    model$obj <- get_objective(tau, rep(lambdaSeq[i], length(tau)), n, m, missInd)
    f_trend <- solve_model(model, solver=solver, y=y)
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
  
  model$obj <- get_objective(tau, lambda, n, m, missInd)
  theta <- solve_model(model, solver=solver, y=y)

  return(list(theta = theta,
              lambda = lambda, 
              BIC = BIC, 
              df = df))
}