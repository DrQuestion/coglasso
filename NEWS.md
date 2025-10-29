# coglasso (development version)

---
editor_options: 
  markdown: 
    wrap: 72
---

# coglasso 1.1.0

## New functions

-   New `bs()` is a wrapping function that single-handedly *builds* the
    multi-omics networks with `coglasso()` and *selects* the best one
    according to the preferred model selection method with
    `select_coglasso()` in a single function call.

-   New `get_network()` is extracts a network in the `igraph` format
    from an object either of class `coglasso` or of class
    `select_coglasso`.

-   New `get_pcor()` extracts a matrix of partial correlations from an
    object either of class `coglasso` or of class `select_coglasso`.

-   New `plot()` can now plot both `coglasso` and `select_coglasso`
    objects. The plots will have color coded nodes and weighted edges.

-   New `select_coglasso()` is a wrapping function to handle all
    possible present (and future) model selection methods. For the
    moment it allows to perform model selection with either *eXtended StARS*, 
    *eXtended Efficient StARS* (see below), or *eBIC*.

-   New `xestars()`, performs *eXtended Efficient StARS*, a
    significantly faster version of *XStARS*.\
    **How much faster?**\
    In our tests, `xestars()` runs 80-90% faster than `xstars()`, even
    more in specific instances.\
    **What features make `xestars()` faster?**\
    First of all, the check for stability that in `xstars()` is
    performed after iterating throughout all the penalty parameters,
    here is implemented as a stopping criterion. Hence, less penalty
    parameters are explored, moreover usually are excluded those that
    lead to denser network (and so to longer network estimations).\
    Second, the use of vectors instead of matrixes to keep track of the
    network variabilities makes the algorithm proceed faster, for the
    former are easier and lighter objects to deal with.\
    Third, a new sampling strategy allows a the computation of as many
    correlation matrixes (the input to `coglasso()`), as the number of
    repetitions of the algorithm only once at the beginning of the
    algorithm. The original strategy performs this every time the
    algorithm switches from the selection of lambda_w to that of a
    lambda_b (which can happen several times). Especially for larger
    data sets, this consists a huge difference.\
    **How do `xstars()` and `xestars()` differ in results?**\
    The impressive increase in speed comes with some minor costs.\
    The different sampling strategy that guarantees not only a faster,
    but also a fairer parameter selection, may lead to different
    selected hyperparameters between the older and the new methodology.

-   New `xstars()` implements the XStARS algorithm seen in the original
    manuscript of *collaborative graphical lasso*. It performs
    stability-based selection of the `c` hyperparameter simultaneously
    with `lambda_w` and `lambda_b`. It substitutes the more primitive
    `stars_coglasso()`, now under deprecation.

## New features and upgrades

-   A **new version** of the *collaborative graphical lasso* algorithm,
    is now able to accept **more than two omics layers**. This new
    version, called *general \|D\|* version, provides the same results
    for two omics layers, but it is slightly slower, so the *general
    \|D\|* algorithm will only be used when necessary. The current
    version has convergence issues for most values of `c`. Hopefully
    this will be fixed by the 2.0.0 release.

-   Added a logo to the package.

-   In `bs()` and `coglasso()`, the generation procedure of`lambda_w`
    and `lambda_b` is now different: the maximum values will be,
    respectively, the highest *within* Pearson's correlation value and
    the highest *between* Pearson's correlation value. Moreover, in
    previous versions the granularity of the search grid increased as
    the values of `lambda_w` and `lambda_b` decreased. As our major
    interest lies in sparser network, this granularity has now been
    inverted.

-   `coglasso()` now outputs an object of `S3` class `coglasso`, while
    all functions whose returned object concerns a selected network, like
    `bs()`,`select_coglasso()`, and all the other selecting functions
    output a `select_coglasso`. Both these classes have related
    `print()` and `plot()` methods.

-   `coglasso()` gains a `lock_lambdas` argument to simulate the
    single penalty parameter-behavior of the original *glasso*. It is
    currently chiefly for testing purposes, so we have not implemented
    any selection procedure for it yet.

## Bug fixes

-   `xstars()` now properly selects `lambda_b`. The selection process
    was never really happening, and we were selecting `lambda_w` twice,
    instead. This will lead to inevitable backward incompatibilities, at
    least of the results of the previous version of `xstars()`, that
    will not be reproducible.

## Deprecations

-   In `coglasso()` and `bs()`, `pX` is being deprecated, will be
    unusable from version 1.2.0 (or 2.0.0). It is now substituted by the
    argument `p`. `p` can take a vector with the dimensions of multiple
    omics layer, as now the package accepts more than two omics layers.

-   `stars_coglasso()` is being deprecated, will be unusable from
    version 1.2.0 (or 2.0.0). Substituted by `xstars()`.

# coglasso 1.0.2

-   Description field of DESCRIPTION now complies the CRAN reviewer's
    comments.

# coglasso 1.0.1

-   Made the Description field of the DESCRIPTION file more explanatory
    and added the official reference to the main method in it.

-   Added the official reference to the main method in the README file.

-   Changed \dontrun{} to \donttest{} in man/stars_coglasso.Rd.

# coglasso 1.0.0

-   Initial CRAN submission.
