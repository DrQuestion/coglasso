test_that("ebic selection works", {
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  sel_cg_ebic <- select_coglasso(cg, method = "ebic", verbose = FALSE)
  expect_equal(c(sel_cg_ebic$sel_index_lw, sel_cg_ebic$sel_index_lb, sel_cg_ebic$sel_index_c), c(3,1,2))
})

test_that("Verbose mode works", {
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_output(select_coglasso(cg, method = "ebic"))
})

test_that("xstars call works", {
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_no_error(select_coglasso(cg, rep_num = 3, verbose = FALSE))
})
