#' Simulate data with peaks and baseline
#' 
#' @param n length of series to generate
#' @export
generate_peaks <- function(n){
  x <- seq(0.5, n, 1)/n
  df <- rpois(1, n/100) + 1
  splineBasis <- ns(x, df=df)
  theta <- rexp(ncol(splineBasis), 1)
  baseline <- splineBasis%*%theta
  
  numberOfPeaks <- rbinom(1, n, .005)
  peakCenters <- round(runif(numberOfPeaks)*(n-2) + 1)
  peakArea <- rnorm(numberOfPeaks, 20, 4)
  peakWidths <- runif(numberOfPeaks)*10 + 2
  
  peaks <- rep(0, length(baseline))
  if (numberOfPeaks > 0){
    for (i in 1:numberOfPeaks){
      peaks <- peaks + 
        peakArea[i]*dnorm(seq(1,n,1), mean=peakCenters[i], sd = peakWidths[i])
    }
  }
  
  noise <- 0.25*rnorm(n)
  y <- peaks + baseline + noise
  df <- data.frame(y=y, x=x, baseline=baseline, peaks = peaks)
  return(df)
}