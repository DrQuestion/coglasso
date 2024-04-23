test_that("bs works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  sel_cg <- bs(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, rep_num = 3, verbose = FALSE)
  expect_equal(c(sel_cg$sel_index_c, sel_cg$sel_index_lw, sel_cg$sel_index_lb), c(1,1,1))
})