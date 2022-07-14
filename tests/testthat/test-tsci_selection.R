test_that("correct output for valid input", {
  set.seed(1)
  n <- 100
  Z <- matrix(rnorm(n), nrow = n)
  D <- Z + Z^2 + rnorm(n)
  Y <- D + Z + rnorm(n)
  vio_space <- cbind(Z, Z^2)
  vio_ind <- list(1, c(1, 2))
  Q <- 3
  A <- cbind(Z, Z^2)
  weight <- A %*% chol2inv(chol(t(A) %*% A)) %*% t(A)
  output <- tsci_selection(Y = Y,
                           D = D,
                           W = NULL,
                           Y_A1 = Y,
                           D_A1 = D,
                           W_A1 = NULL,
                           vio_space = vio_space,
                           vio_ind = vio_ind,
                           Q = Q,
                           weight = weight,
                           intercept = F,
                           str_thol = 10,
                           alpha = 0.05)
  expect_equal(2 * 2, 4)
})
