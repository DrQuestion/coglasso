---
editor_options: 
  markdown: 
    wrap: 72
---

# coglasso (development version)

-   Began the deprecation process of `stars_coglasso()`, renamed to the
    more accurate `xstars()`. Made also some implementation improvements
    to `xstars()` that made it faster.

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
    significantly faster and more memory-efficient version of *XStARS*.\
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
    **What features make `xestars()` more memory-efficient?\
    **`xestars()` does not keep in memory for longer than strictly
    necessary the "merged" matrixes that store the average stability of
    each edge, without accumulating them as in the original `xstars()`.
    This, unless not explicitely requested by setting the parameter
    `light` to `FALSE`.\
    **How do `xstars()` and `xestars()` differ in results?\
    **The impressive increase in speed and more memory-efficiency comes
    with some minor costs.\
    First of all, some objects returned by `xstars()` are not returned
    if not explecetely requested, and in general they come in a smaller
    amount (see question above).\
    Second, the different sampling strategy may guarantee not only a
    faster, but also a fairer parameter selection, as they are all
    selected from the same subsamplings. This may lead to different
    selected hyperparameters between the older and the new methodology.

-   Created the classes `coglasso` and `select_coglasso`, with related
    `print()` and `plot()` methods. These former is returned when using
    `coglasso()`, while the latter when using any function that performs
    model selection as its final step (*i.e.*
    `bs()`,`select_coglasso()`, and all the model selection functions).

-   Created a plotting module. Now every function of the package that
    produces or selects a network generates an object that can be
    directly be plotted with `plot()`.

-   Implemented a **new version** of the *collaborative graphical lasso*
    algorithm, now able to accept **more than two omics layers**. These
    new version, called *general \|D\|* version, is perfectly compatible
    with and as time consuming as the previous one when only two layers
    are given as input.

-   Began the deprecation process of the argument `pX` of the functions
    `coglasso()` and `bs()`, substituted by the argument `p`, that can
    take a vector with the dimensions of multiple omics layer, as now
    the package accepts more than two omics layers.

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
