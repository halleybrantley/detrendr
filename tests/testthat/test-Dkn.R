context("Discrete Derivatives")

test_that("Dk(n) produces expected output",{
  n <- 10
  k <- 4
  expect_that(get_Dk(n, k), is_equivalent_to(get_Dk_R(n,k)))
})
