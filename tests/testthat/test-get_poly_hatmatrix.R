test_that("correct output for valid input", {
  set.seed(1)
  n <- 100
  p <- 2
  df <- readRDS("df_get_poly_hatmatrix.rds")
  Z <- matrix(df$Z[seq_len(n * p)], nrow = n, ncol = p)
  df <- cbind(1, Z)
  params_list <- list(3, 2)
  A <- cbind(1, Z[, 1], Z[, 1]^2, Z[, 1]^3, Z[, 2], Z[, 2]^2)
  output <- get_poly_hatmatrix(df_treatment_A1 = df, params = params_list)
  expect_equal(output$weight,
               A %*% chol2inv(chol(t(A) %*% A)) %*% t(A))

})
