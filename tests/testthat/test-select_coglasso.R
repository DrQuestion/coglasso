test_that("ebic selection works", {
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  sel_cg_ebic <- select_coglasso(cg, method = "ebic", verbose = FALSE)
  expect_equal(c(sel_cg_ebic$sel_index_lw, sel_cg_ebic$sel_index_lb, sel_cg_ebic$sel_index_c), c(3,1,2))
})

test_that("Verbose mode works", {
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_output(select_coglasso(cg, method = "ebic"))
})

test_that("xestars call works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_no_error(select_coglasso(cg, rep_num = 3, verbose = FALSE))
})

test_that("xstars call works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_no_error(select_coglasso(cg, method = "xstars", rep_num = 3, verbose = FALSE))
})

test_that("print.select_coglasso works with all functions", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  expect_snapshot({
    cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
    sel_cg <- select_coglasso(cg, method = "ebic", verbose = FALSE)
    print(sel_cg)
    sel_cg <- select_coglasso(cg, method = "xestars", rep_num = 3, verbose = FALSE)
    print(sel_cg)
    sel_cg <- xstars(cg, rep_num = 3, verbose = FALSE)
    print(sel_cg)
    sel_cg <- bs(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, rep_num = 3, verbose = FALSE)
    print(sel_cg)
  })
  
})