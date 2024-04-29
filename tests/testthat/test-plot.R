test_that("plot.select_coglasso and plot.coglasso work", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  # plot.coglasso
  cg <- coglasso(multi_omics_sd_small, pX = 14, nlambda_w = 15, nlambda_b = 15,
               nc = 2, lambda_w_min_ratio = 0.6, verbose = FALSE)
  expect_no_error(plot(cg, index_c = 1, index_lw = 6, index_lb = 6))
  expect_no_error(plot(cg, index_c = 1, index_lw = 6, index_lb = 6, node_labels = FALSE, hide_isolated = FALSE))
  # plot.select_coglasso
  sel_cg <- select_coglasso(cg, rep_num = 4, verbose = FALSE)
  expect_no_error(plot(sel_cg))
  expect_no_error(plot(sel_cg, node_labels = FALSE, hide_isolated = FALSE))
})
