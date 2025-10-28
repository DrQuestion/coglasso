test_that("warning raised when not implemented method is selected", {
  withr::with_seed(42, {
    cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 2, nlambda_b = 2, nc = 1, verbose = FALSE)
    expect_warning(select_coglasso(cg, method = "aic", rep_num = 3, verbose = FALSE))
  })
})

test_that("ebic selection works", {
  cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  sel_cg_ebic <- select_coglasso(cg, method = "ebic", verbose = FALSE)
  expect_equal(c(sel_cg_ebic$sel_index_lw, sel_cg_ebic$sel_index_lb, sel_cg_ebic$sel_index_c), c(3,1,1))
})

test_that("Verbose mode works", {
  cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_output(select_coglasso(cg, method = "ebic"))
})

test_that("xestars call works", {
  withr::with_seed(42, {
    cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
    expect_no_error(select_coglasso(cg, rep_num = 3, verbose = FALSE))
  })
})

test_that("xstars call works", {
  withr::with_seed(42, {
    cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
    expect_no_error(select_coglasso(cg, method = "xstars", rep_num = 3, verbose = FALSE))
  })
})

test_that("print.select_coglasso works with all functions", {
  withr::with_seed(42, {
    expect_snapshot({
      cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
      sel_cg <- select_coglasso(cg, method = "ebic", verbose = FALSE)
      print(sel_cg)
      sel_cg <- select_coglasso(cg, method = "xestars", rep_num = 3, verbose = FALSE)
      print(sel_cg)
      sel_cg <- xstars(cg, rep_num = 3, verbose = FALSE)
      print(sel_cg)
      sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, rep_num = 3, verbose = FALSE)
      print(sel_cg)
    })
  })
})