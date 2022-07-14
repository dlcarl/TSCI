test_that("correct output for valid input", {
  set.seed(1)
  n <- 100
  p <- 2
  Z <- matrix(rnorm(n * p), nrow = n)
  vio_space <- list(Z, Z^2, Z^3)
  list_vio_space <- build_vio_space_candidates(vio_space, create_nested_sequence = TRUE)
  expect_equal(list_vio_space$vio_space, cbind(Z, Z^2, Z^3))
  expect_equal(lapply(list_vio_space$vio_ind, FUN = function(x) sort(x)), list(c(1:2), c(1:4), c(1:6)))
  expect_equal(list_vio_space$Q, 4)
  expect_equal(list_vio_space$nested_sequence, TRUE)

  list_vio_space <- build_vio_space_candidates(vio_space, create_nested_sequence = FALSE)
  expect_equal(list_vio_space$vio_space, cbind(Z, Z^2, Z^3))
  expect_equal(lapply(list_vio_space$vio_ind, FUN = function(x) sort(x)), list(c(1:2), c(3:4), c(5:6)))
  expect_equal(list_vio_space$Q, 4)
  expect_equal(list_vio_space$nested_sequence, FALSE)
})
