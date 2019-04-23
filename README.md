
<!-- README.md is generated from README.Rmd. Please edit that file -->
detrendr
========

The goal of detrendr is to estimate smooth quantile trends in evenly spaced series. Missing values are allowed.

Installation
============

You can install the development version of detrendr from [Github](https://github.com/halleybrantley/detrendr) with:

``` r
library(devtools)
install_github("halleybrantley/detrendr")
#> Skipping install of 'detrendr' from a github remote, the SHA1 (c277bd60) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

Example
-------

Estimate the 5th and 20th quantile trends using a fixed value of the smoothing parameter *λ*.

``` r
library(detrendr)
#> Loading required package: Matrix
n <- 100
x <- seq(1, n, 1)
y <- sin(x*2*pi/n) + rnorm(n, 0, .4)
lambda <- 10
k <- 3
tau <- c(0.05, .2)
trend <- get_trend(y, tau, lambda, k)
#> Using same lambda for all quantiles
plot(y~x, type="l", col="grey")
lines(trend[,1]~x, col="red")
lines(trend[,2]~x, col="blue")
```

<img src="man/figures/README-example1-1.png" width="70%" />

Use eBIC criterion to choose smoothing parameter:

``` r

trend_fit <- get_trend_BIC(y, tau, k, plot_lambda = TRUE)
#> Using same lambda for all quantiles
```

<img src="man/figures/README-BIC-1.png" width="70%" />

``` r
trend <- trend_fit$trend
plot(y~x, type="l", col="grey")
lines(trend[,1]~x, col="red")
lines(trend[,2]~x, col="blue")
```

<img src="man/figures/README-BIC-2.png" width="70%" />
