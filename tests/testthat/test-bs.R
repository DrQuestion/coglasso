test_that("bs works", {
  withr::with_seed(42, {
    sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, rep_num = 3, verbose = FALSE)
  })
  expect_equal(c(sel_cg$sel_index_c, sel_cg$sel_index_lw, sel_cg$sel_index_lb), c(1,1,1))
})

test_that("whole bs/xstars run with 1 combination of hpars works", {
  withr::with_seed(42, {
    sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 1, nlambda_b = 1, nc = 1, rep_num = 3, method = "xstars", verbose = FALSE)
  })
  expect_equal(c(sel_cg$sel_index_c, sel_cg$sel_index_lw, sel_cg$sel_index_lb), c(1,1,1))
})

test_that("whole bs/xestars run with 1 combination of hpars works", {
  withr::with_seed(42, {
    sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 1, nlambda_b = 1, nc = 1, rep_num = 3, method = "xestars", verbose = FALSE)
  })
  expect_equal(c(sel_cg$sel_index_c, sel_cg$sel_index_lw, sel_cg$sel_index_lb), c(1,1,1))
})
