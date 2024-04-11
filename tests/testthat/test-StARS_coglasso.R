test_that("stars_coglasso works", {
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_no_error(stars_coglasso(cg, rep_num = 3, verbose = FALSE))
})

test_that("Verbose mode works", {
  cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_output(stars_coglasso(cg, rep_num = 3))
})
