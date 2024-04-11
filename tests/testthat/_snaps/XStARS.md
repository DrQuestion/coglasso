# stars_coglasso is deprecated

    Code
      cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 2,
        verbose = FALSE)
      expect_no_error(stars_coglasso(cg, rep_num = 3, verbose = FALSE))
    Condition
      Warning:
      `stars_coglasso()` was deprecated in coglasso 1.1.0.
      i Please use `xstars()` instead.

