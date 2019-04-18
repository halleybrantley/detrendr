
# Criteria for evaluating the smoothness parameter
#'
#' \code{get_criteria}
#'
#' @param criteria label of criteria to be used, must be one of "eBIC", "SIC", "valid"
#' @param f_trend matrix of fitted quantile trend(s)
#' @param y observed data vector
#' @param tau quantile vector
#' @param D discrete differencing matrix (for SIC and eBIC)
#' @param df_tol tolerance for determining degrees of freedom (Dtheta > df_tol)
#' @param gamma parameter for eBIC
#' @param validID index of data to be used for validation (for valid method) 
#' @param yValid validation data (for valid method) 
#' @export
get_criteria <- function(criteria, f_trend, y, tau, 
                         D = NULL, df_tol=1e-9, gamma=1, 
                         validID = NULL, yValid = NULL){
  
  if (criteria == "eBIC" || criteria == "SIC") {
    n <- length(y)
    missInd <- which(is.na(y))
    if (sum(missInd) > 0){
      resid_trend <- checkloss(y[-missInd]-f_trend[-missInd,,drop=FALSE], tau)
    } else {
      resid_trend <- checkloss(y-f_trend, tau)
    }
    df <- Matrix::colSums(abs(D%*%f_trend) > df_tol) 
    if (criteria == "eBIC"){
      scale_param <- 0.5 - abs(0.5-tau)
      BIC <- 2*colSums(resid_trend, na.rm=T)/scale_param + log(n)*df +
        2*gamma*log(choose(nrow(D), df))
      if (any(sapply(BIC, is.infinite))) {
        H <- function(x) {
          return(x*log(1/x) + (1-x)*log(1/(1-x)))
        }
        BIC <- 2*colSums(resid_trend, na.rm=T)/scale_param + log(n)*df +
          2*gamma*nrow(D)*H(df/nrow(D))
      }
    } else {
      BIC <- log(colMeans(resid_trend, na.rm=T)) + log(n)*df/(2*n)
    }
  } else if (criteria == "valid"){
    missInd <- is.na(yValid)
    if (sum(missInd) > 0){
      BIC <- colMeans(checkloss(
        yValid[-missInd]-f_trend[validID[-missInd],,drop=FALSE], tau))
    } else {
      BIC <- colMeans(checkloss(
        yValid-f_trend[validID,,drop=FALSE], tau))
    }
    
    df <- NA
  } 
  
  return(list(BIC = BIC, df=df))
} 