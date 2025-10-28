test_that("get_network stop conditions work", {
  a <- 0
  class(a) <- "coglasso"
  expect_error(get_network(a), 'input object is of class "coglasso"')
  a <- list()
  class(a) <- "select_coglasso"
  a$method <- "ebic"
  expect_error(get_network(a, index_c = 1), "the others must be specified")
  a$method <- "xstars"
  expect_error(get_network(a, index_c = 1), "the others must be specified")
})

test_that("get_network from ebic works", {
  sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, method = "ebic", verbose = FALSE)
  expect_s3_class(get_network(sel_cg), "igraph")
  expect_s3_class(get_network(sel_cg, index_c = 2, index_lw = 2, index_lb = 2), "igraph")
})

test_that("get_network from coglasso works", {
  cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_s3_class(get_network(cg, index_c = 2, index_lw = 2, index_lb = 2), "igraph")
})

test_that("get_network from xstars/xestars works", {
  withr::with_seed(42, {
    sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, rep_num = 3, verbose = FALSE)
    expect_s3_class(get_network(sel_cg), "igraph")
    expect_s3_class(get_network(sel_cg, index_c = 2, index_lw = 2, index_lb = 2), "igraph")
  })
})

test_that("get_network without colnames works", {
  a <- multi_omics_sd_micro
  colnames(a) <- NULL
  cg <- coglasso(a, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_s3_class(get_network(cg, index_c = 2, index_lw = 2, index_lb = 2), "igraph")
})

test_that("get_pcor stop conditions work", {
  a <- 0
  class(a) <- "coglasso"
  expect_error(get_pcor(a), 'input object is of class "coglasso"')
  a <- list()
  class(a) <- "select_coglasso"
  a$method <- "ebic"
  expect_error(get_pcor(a, index_c = 1), "the others must be specified")
  a$method <- "xstars"
  expect_error(get_pcor(a, index_c = 1), "the others must be specified")
})

test_that("get_pcor from ebic works", {
  sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, method = "ebic", verbose = FALSE)
  expect_equal(dim(get_pcor(sel_cg)), c(6,6))
  expect_equal(dim(get_pcor(sel_cg, index_c = 2, index_lw = 2, index_lb = 2)), c(6,6))
})

test_that("get_pcor from coglasso works", {
  cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_equal(dim(get_pcor(cg, index_c = 2, index_lw = 2, index_lb = 2)), c(6,6))
})

test_that("get_pcor from xstars/xestars works", {
  withr::with_seed(42, {
    sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, rep_num = 3, verbose = FALSE)
    expect_equal(dim(get_pcor(sel_cg)), c(6,6))
    expect_equal(dim(get_pcor(sel_cg, index_c = 2, index_lw = 2, index_lb = 2)), c(6,6))
  })
})

test_that("get_pcor without colnames works", {
  a <- multi_omics_sd_micro
  colnames(a) <- NULL
  cg <- coglasso(a, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2, verbose = FALSE)
  expect_equal(dim(get_pcor(cg, index_c = 2, index_lw = 2, index_lb = 2)), c(6,6))
})
