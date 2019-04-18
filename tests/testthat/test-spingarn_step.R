context("Spingarn Step")



test_that("Multiple Spingarn steps produce expected result", {
  require(Matrix)
  set.seed(12345)
  n <- 1e2
  x <- seq(1/n, 1, length.out=n)
  f <- 2*(x + 2)^2 + 3*cos(3*pi*x)
  tau <- 1e4
  g <-100*exp(-tau*(x-0.5)^2)
  y <- f + g + rnorm(n)
  k <- 3
  D <- get_Dk(n, k)
  M <- diag(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)
  lambda <- 1
  tau <- 0.01
  step <- 1
  theta <- y
  eta <- as.numeric(matrix(D %*% theta))
  multi_step <- spingarn_multi_step(theta, eta, y, D, cholM, lambda, tau, 
                                    step, 5, k)
  
  expect_that(spingarn_multi_step_R(theta, eta, y, D, M, lambda, tau, step,5), 
              is_equivalent_to(
                spingarn_multi_step(theta, eta, y, D, cholM, lambda, tau, 
                                    step, 5, k)))
  expect_that(theta, equals(y))
})
  
