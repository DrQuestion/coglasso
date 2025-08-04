test_that("plot.select_coglasso and plot.coglasso work", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  # plot.coglasso
  cg <- coglasso(multi_omics_sd_small, p = 14, nlambda_w = 5, nlambda_b = 5,
               nc = 2, lambda_w_min_ratio = 0.4, lambda_b_min_ratio = 0.5, verbose = FALSE)
  expect_no_error(plot(cg, index_c = 2, index_lw = 1, index_lb = 5))
  expect_no_error(plot(cg, index_c = 2, index_lw = 1, index_lb = 5, node_labels = FALSE, hide_isolated = FALSE))
  # plot.select_coglasso
  sel_cg <- select_coglasso(cg, rep_num = 3, verbose = FALSE)
  expect_no_error(plot(sel_cg))
  expect_no_error(plot(sel_cg, node_labels = FALSE, hide_isolated = FALSE))
})

test_that("get_fillsframes works with more than 6 data sets", {
  expect_no_error(get_fillsframes(7))
})

