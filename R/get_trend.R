#' Get quantile trends
#'
#' \code{get_trend} Returns the quantile trend matrix
#'
#' @param y observed data, should be equally spaced, may contain NA
#' @param tau quantile levels at which to evaluate trend
#' @param lambda penalty paramter controlling smoothness
#' @param k order of differencing
#' @importFrom utils installed.packages
#' @examples 
#' n <- 100
#' x <- seq(1, n, 1)
#' y <- sin(x*2*pi/n) + rnorm(n, 0, .4)
#' lambda <- 10
#' k <- 3
#' tau <- c(0.05, .2)
#' trend <- get_trend(y, tau, lambda, k)
#' plot(y~x, type = "l")
#' lines(trend[,1]~x, col = "red")
#' @export
get_trend <- function(y, tau, lambda, k){
  mean_y <- mean(y, na.rm=T)
  sd_y <- stats::sd(y, na.rm=T)
  y <- as.numeric(scale(y))*200
  model <- get_model(y, tau, lambda, k)
  pkgs <- installed.packages()[,"Package"]
  if("gurobi" %in% pkgs){
    solver <- "gurobi"
  } else if ("Rglpk" %in% pkgs){
    solver <- "Rglpk"
  } else {
    solver <- "lpSolve"
  }
  trend <- solve_model(model, solver, y)
  return(trend/200*sd_y + mean_y)
}

