context("Proximal Mapping")

test_that("prox_quantile produces same output as R function",{
  n <- 1e3
  w <- seq(-3, 3, length.out=n)
  tau <- 0.5
  alpha <- 2
  
  expect_that(prox_quantile(w, tau, alpha), 
              is_equivalent_to(prox_quantile_R(w, tau, alpha)))
})

test_that("prox_f1 produces same output as R function",{
  set.seed(1234)
  theta <- rnorm(100)
  y <- theta + rnorm(100, .01)
  expect_that(prox_f1(theta, y, 0.05, 1), 
              is_equivalent_to(prox_f1_R(theta, y, 0.05, 1)))
})

test_that("prox_f2 produces same output as R function", {
  eta <- seq(-3, 3, length.out = 1000)
  lambda <- 1
  step <- 1
  expect_that(prox_f2(eta, lambda, step), 
              is_equivalent_to(prox_f2_R(eta, lambda, step)))
  expect_that(prox_f2(eta, lambda), 
              is_equivalent_to(prox_f2_R(eta, lambda)))
})

test_that("prox produces same output as R function", {
  set.seed(1234)
  theta <- rnorm(100)
  y <- theta + rnorm(100, .01)
  eta <- seq(-3, 3, length.out = 1000)
  lambda <- 1
  step <- 1
  tau <- 0.05
  expect_that(prox_test(theta, eta,  y, lambda, tau, step), 
              is_equivalent_to(prox_R(theta, eta,  y, lambda, tau, step)))
  expect_that(prox_test(theta, eta,  y, lambda), 
              is_equivalent_to(prox_R(theta, eta,  y, lambda)))
})





