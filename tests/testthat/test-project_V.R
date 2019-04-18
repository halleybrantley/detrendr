context("Projection onto V")

test_that("Projection onto V produces expected result", {
  require(Matrix)
  set.seed(12345)
  k <- 3
  n <- 1e2
  D <- get_Dk(n,k)
  M <- diag(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)
  theta <- rnorm(n)
  eta <- as.numeric(D%*%theta) + 0.01*rnorm(n-k)
  proj1 <- project_V_R(theta, eta, D, M)
  theta1 <- as.numeric(proj1[[1]])
  eta1 <- as.numeric(proj1[[2]])
  # Alters theta and eta in place
  project_V(theta, eta, D, cholM, k)
  expect_that(theta1, 
             is_equivalent_to(theta))
  expect_that(eta1, 
              is_equivalent_to(eta))
  
})

