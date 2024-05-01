test_that("xstars works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  sel_cg <- xstars(cg, rep_num = 3, verbose = FALSE)
  expect_equal(c(sel_cg$sel_index_c, sel_cg$sel_index_lw, sel_cg$sel_index_lb), c(2,1,1))
})

test_that("Verbose mode of xstars works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_output(xstars(cg, rep_num = 3))
})

test_that("xestars works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  sel_cg <- xestars(cg, rep_num = 3, verbose = FALSE)
  expect_equal(c(sel_cg$sel_index_c, sel_cg$sel_index_lw, sel_cg$sel_index_lb), c(1,1,1))
  sel_cg <- xestars(cg, rep_num = 3, light = FALSE, old_sampling = TRUE, verbose = FALSE)
  expect_equal(c(sel_cg$sel_index_c, sel_cg$sel_index_lw, sel_cg$sel_index_lb), c(2,1,1))
})

test_that("Verbose mode of xestars works", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_output(xestars(cg, rep_num = 3))
})

test_that("stars_coglasso is deprecated", {
  old_seed <- get0(".Random.seed", envir = .GlobalEnv)
  on.exit({
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  set.seed(42)
  expect_snapshot({
    cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
    expect_no_error(stars_coglasso(cg, rep_num = 3, verbose = FALSE))
  })
})
