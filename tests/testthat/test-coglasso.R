test_that("coglasso works", {
  expect_no_error(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3, verbose = FALSE))
})

test_that("nlambda_w works", {
  n <- 3
  expect_length(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, verbose = FALSE)$lambda_w, n)
})

test_that("nlambda_b works", {
  n <- 3
  expect_length(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, verbose = FALSE)$lambda_b, n)
})

test_that("nc works", {
  n <- 3
  expect_length(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, verbose = FALSE)$c, n)
})

test_that("lambda_w_max works", {
  n <- 3
  lambda_w_max <- 0.8
  expect_equal(max(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, lambda_w_max = lambda_w_max, verbose = FALSE)$lambda_w), lambda_w_max)
})

test_that("lambda_b_max works", {
  n <- 3
  lambda_b_max <- 0.8
  expect_equal(max(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, lambda_b_max = lambda_b_max, verbose = FALSE)$lambda_b), lambda_b_max)
})

test_that("c_max works", {
  n <- 3
  c_max <- 20
  expect_equal(max(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, c_max = c_max, verbose = FALSE)$c), c_max)
})

test_that("lambda_w_min_ratio works", {
  n <- 3
  lambda_w_min_ratio <- 0.2
  actual_min_lambda_w <- min(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, lambda_w_min_ratio = lambda_w_min_ratio, verbose = FALSE)$lambda_w)
  actual_max_lambda_w <- max(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, lambda_w_min_ratio = lambda_w_min_ratio, verbose = FALSE)$lambda_w)
  expect_equal(actual_min_lambda_w/actual_max_lambda_w, lambda_w_min_ratio)
})

test_that("lambda_b_min_ratio works", {
  n <- 3
  lambda_b_min_ratio <- 0.2
  actual_min_lambda_b <- min(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, lambda_b_min_ratio = lambda_b_min_ratio, verbose = FALSE)$lambda_b)
  actual_max_lambda_b <- max(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, lambda_b_min_ratio = lambda_b_min_ratio, verbose = FALSE)$lambda_b)
  expect_equal(actual_min_lambda_b/actual_max_lambda_b, lambda_b_min_ratio)
})

test_that("c_min_ratio works", {
  n <- 3
  c_min_ratio <- 0.2
  actual_min_c <- min(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, c_min_ratio = c_min_ratio, verbose = FALSE)$c)
  actual_max_c <- max(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = n, nlambda_b = n, nc = n, c_min_ratio = c_min_ratio, verbose = FALSE)$c)
  expect_equal(actual_min_c/actual_max_c, c_min_ratio)
})

test_that("single lambda_w works and overrides nlambda_w", {
  n <- 3
  lambda_w <- 0.5
  actual_lambda_w <- coglasso(multi_omics_sd_micro, pX = 4, lambda_w = lambda_w, nlambda_w = n, nlambda_b = n, nc = n, verbose = FALSE)$lambda_w
  expect_length(actual_lambda_w, 1)
})

test_that("single lambda_b works and overrides nlambda_b", {
  n <- 3
  lambda_b <- 0.5
  actual_lambda_b <- coglasso(multi_omics_sd_micro, pX = 4, lambda_b = lambda_b, nlambda_w = n, nlambda_b = n, nc = n, verbose = FALSE)$lambda_b
  expect_length(actual_lambda_b, 1)
})

test_that("single c works and overrides nc", {
  n <- 3
  c <- 0.5
  actual_c <- coglasso(multi_omics_sd_micro, pX = 4, c = c, nlambda_w = n, nlambda_b = n, nc = n, verbose = FALSE)$c
  expect_length(actual_c, 1)
})

test_that("Verbose mode works", {
  expect_output(coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3))
})

test_that("Covariance as input works", {
  S <- stats::cor(multi_omics_sd_micro)
  expect_no_error(coglasso(S, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3, verbose = FALSE))
  expect_output(coglasso(S, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3))
})

test_that("Default run works", {
  expect_no_error(coglasso(multi_omics_sd_micro, pX = 4, verbose = FALSE))
})

test_that("gen_hpar getc parameter works", {
  S <- stats::cor(multi_omics_sd_micro)
  d <- dim(S)[[1]]
  hpars <- gen_hpars(S, d, getc = TRUE)
  expect_equal(colnames(hpars[[1]])[1], "c")
})

test_that("print.coglasso works", {
  expect_snapshot({
    cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
    print(cg)
  })
})