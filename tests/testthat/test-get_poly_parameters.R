test_that("selected polynomials fit model well using 'backfitting'", {
  n <- 100
  p <- 1
  df <- readRDS("df_get_poly_parameters.rds")
  Z <- matrix(df$Z[seq_len(n * p)], nrow = n, ncol = p)
  Y <- 1 + Z + Z^2 + df$error[seq_len(n)]
  df <- data.frame(Y, Z)
  params_list <- list(seq(1, 4, by = 1))
  poly_CV1 <- get_poly_parameters(df_treatment = df,
                                 params_list = params_list,
                                 order_selection_method = "backfitting",
                                 max_iter = 100,
                                 conv_tol = 10^-6,
                                 gcv = TRUE,
                                 nfolds = 5)
  expect_lt(poly_CV1$mse, 1.2)
})

test_that("polynomials are correct for one variable using 'grid search'", {
  n <- 100
  p <- 1
  df <- readRDS("df_get_poly_parameters.rds")
  Z <- matrix(df$Z[seq_len(n * p)], nrow = n, ncol = p)
  Y <- 1 + Z + Z^2 + df$error[seq_len(n)]
  df <- data.frame(Y, Z)
  params_list <- list(seq(1, 4, by = 1))
  poly_CV1 <- get_poly_parameters(df_treatment = df,
                                 params_list = params_list,
                                 order_selection_method = "grid search",
                                 max_iter = 100,
                                 conv_tol = 10^-6,
                                 gcv = TRUE,
                                 nfolds = 5)
  expect_lt(poly_CV1$mse, 1.2)
})

test_that("polynomials are correct for more than one variable using 'backfitting'", {
  n <- 100
  p <- 3
  df <- readRDS("df_get_poly_parameters.rds")
  Z <- matrix(df$Z[seq_len(n * p)], nrow = n, ncol = p)
  Y <- 1 + Z[, 1] + Z[, 1]^2 + Z[, 2]^3 + Z[, 3] + df$error[seq_len(n)]
  df <- data.frame(Y, Z)
  params_list <- list(seq(1, 3, by = 1), seq(1, 3, by = 1), seq(1, 3, by = 1))
  poly_CV1 <- get_poly_parameters(df_treatment = df,
                                  params_list = params_list,
                                  order_selection_method = "backfitting",
                                  max_iter = 100,
                                  conv_tol = 10^-6,
                                  gcv = TRUE,
                                  nfolds = 5)
  expect_lt(poly_CV1$mse, 1.2)
})

test_that("polynomials are correct for more than one variable using 'grid search'", {
  n <- 100
  p <- 3
  df <- readRDS("df_get_poly_parameters.rds")
  Z <- matrix(df$Z[seq_len(n * p)], nrow = n, ncol = p)
  Y <- 1 + Z[, 1] + Z[, 1]^2 + Z[, 2]^3 + Z[, 3] + df$error[seq_len(n)]
  df <- data.frame(Y, Z)
  params_list <- list(seq(1, 3, by = 1), seq(1, 3, by = 1), seq(1, 3, by = 1))
  poly_CV1 <- get_poly_parameters(df_treatment = df,
                                  params_list = params_list,
                                  order_selection_method = "grid search",
                                  max_iter = 100,
                                  conv_tol = 10^-6,
                                  gcv = TRUE,
                                  nfolds = 5)
  expect_lt(poly_CV1$mse, 1.2)
})
