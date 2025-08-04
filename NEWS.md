---
editor_options: 
  markdown: 
    wrap: 72
---

# coglasso (development version)

-   Changed the `xstars()` algorithm to perform the stability-based selection of
    the `c` hyperparameter simultaneously with `lambda_w` and `lambda_b`.

-   Began the deprecation process of `stars_coglasso()`, renamed to the
    more accurate `xstars()`. Made also some implementation improvements
    to `xstars()` that made it faster.

-   Corrected major bug in `xstars()`, as the selection of `lambda_b`
    was never really happening, and we were selecting `lambda_w` twice,
    instead. This will lead to inevitable backward incompatibilities, at
    least of the results of the previous version of `xstars()`, that
    will not be reproducible.

-   Implemented `select_coglasso()` for handling all possible present
    (and future) model selection methods from a single wrapping
    function. For the moment it allows to perform model selection with
    either *eXtended StARS*, *eXtended* *Efficient StARS* (see below),
    or *eBIC*.

-   Implemented `bs()` a single wrapping function to *build* the
    multi-omics networks with `coglasso()` and *select* the best one
    according to the preferred model selection method with
    `select_coglasso()` in a single function call.

-   Implemented `xestars()`, performing *eXtended Efficient StARS*, a
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
    **How do `xstars()` and `xestars()` differ in results?\
    **The impressive increase in speed comes
    with some minor costs.\
    The different sampling strategy that guarantees not only a
    faster, but also a fairer parameter selection, may lead to different
    selected hyperparameters between the older and the new methodology.

-   Created the classes `coglasso` and `select_coglasso`, with related
    `print()` and `plot()` methods. These former is returned when using
    `coglasso()`, while the latter when using any function that performs
    model selection as its final step (*i.e.*
    `bs()`,`select_coglasso()`, and all the model selection functions).

-   Created a plotting module. Now every function of the package that
    produces or selects a network generates an object that can be
    directly be plotted with `plot()`.

-   Began the implementation of a **new version** of the *collaborative
    graphical lasso* algorithm, now able to accept **more than two omics
    layers**. This new version, called *general \|D\|* version, provides
    the same results for two omics layers, but it is slightly slower, so
    the *general \|D\|* algorithm will only be used when necessary. The
    current version has convergence issues for most values of `c`.
    Hopefully this will be fixed before the 2.0 release.

-   Began the deprecation process of the argument `pX` of the functions
    `coglasso()` and `bs()`, substituted by the argument `p`, that can
    take a vector with the dimensions of multiple omics layer, as now
    the package accepts more than two omics layers.

-   Added a logo to the package, but it will be changed to a white
    background.

-   Smaller functional bugs have been corrected.

-   A new parameter for `coglasso()`, `lock_lambdas`, has been
    introduced to simulate the single penalty parameter-beahavior of the
    original *glasso*. It is currently chiefly for testing purposes, so
    we have not implemented any selection procedure for it, yet.

-   Changed the generation procedures for `lambda_w` and `lambda_b`: the
    maximum values will be, respectively, the highest *within* Pearson's
    correlation value and the highest *between* Pearson's correlation
    value. Moreover, in previous versions the granularity of the search 
    grid increased as the values of `lambda_w` and `lambda_b` decreased. 
    As our major interest lies in sparser network, this granularity has 
    now been inverted. 
    
-   Implemented the new functions `get_network()` and `get_pcor()`. These
    two functions extract, respectively, a network in the `igraph` format
    or a matrix of partial correlations from an object either of class
    `coglasso` or of class `select_coglasso`.

# coglasso 1.0.2

-   Reformatted the Description field of DESCRIPTION file according to
    CRAN reviewer's comments.

# coglasso 1.0.1

-   Made the Description field of the DESCRIPTION file more explanatory
    and added the official reference to the main method in it.

-   Added the official reference to the main method in the README file.

-   Changed \dontrun{} to \donttest{} in man/stars_coglasso.Rd.

# coglasso 1.0.0

-   Initial CRAN submission.
