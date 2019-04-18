#' Get evaluation metrics for classification
#'
#' \code{get_metric} Returns classification metric requested
#'
#' @param trend estimated trend vector
#' @param y vector of observations
#' @param signal vector of true signal classification
#' @param threshold value for threshold defining signal
#' @param metric choice of "missclass", "recall", "CAA", "precision" or "F1"
#' @export
get_metric <- function(trend, y, signal, threshold, metric){
  y_adj <- y - trend
  signal_hat <- as.numeric(y_adj > threshold)
  if (metric == "missclass"){
    return(mean(abs(signal-signal_hat), na.rm=T))
  } else {
    true_pos <- sum(signal_hat==1 & signal==1, na.rm=T)
    recall <- true_pos/sum(signal == 1, na.rm=T)
    if (metric == "recall"){
      return(recall)
    } else if (metric == "CAA"){
      recall2 <- sum(signal_hat==0 & signal==0, na.rm=T) / 
        sum(signal == 0, na.rm=T)
      caa <- (recall + recall2)/2
      return(caa)
    }
    precision <- true_pos/sum(signal_hat == 1, na.rm=T)
    if (metric == "precision"){
      return(precision)
    } else if (metric == "F1"){
      F1 <- 2 * (precision*recall)/(precision + recall)
      return(F1)
    } else {
      stop("Metric must be in missclass, recall, precistion, F1, CAA")
    }
  }
}