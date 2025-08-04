#' Select the best `coglasso` network
#'
#' `select_coglasso()` selects the best combination of hyperparameters given to
#' `coglasso()` according to the selected model selection method. The three
#' availble options that can be set for the argument `method` are  "xstars", 
#' "xestars" and "ebic".
#'
#' `select_coglasso()` provides three model selection strategies:
#' \itemize{
#'   \item "xstars" uses *eXtended StARS* (*XStARS*) selecting the most stable, yet sparse
#'     network. Stability is computed upon network estimation from multiple subsamples of the
#'     multi-omics data set, allowing repetition. Subsamples are collected for a
#'     fixed amount of times (`rep_num`), and with a fixed proportion of the total
#'     number of samples (`stars_subsample_ratio`). See [xstars()] for more
#'     information on the methodology.
#'   \item "xestars" uses *eXtended Efficient StARS* (*XEStARS*), a significantly
#'     faster version of *XStARS*. It could produce marginally different results
#'     to "xstars" due to a different sampling strategy. See [xestars()] for 
#'     more information on the methodology.
#'   \item "ebic" uses the *extended Bayesian Information*
#'     *Criterion* (*eBIC*) selecting the network that minimizes it. `gamma` sets the
#'     wait given to the extended component, turning the model selection method to
#'     the standard *BIC* if set to 0.
#' }
#' 
#' @inherit xestars
#' @param method The model selection method to select the best combination of
#'   hyperparameters. The available options are "xstars", "xestars" and "eBIC". 
#'   Defaults to "xestars".
#' @param ebic_gamma The \eqn{\gamma} tuning parameter for *eBIC* selection, to set
#'   between 0 and 1. When set to 0 one has the standard *BIC*. Defaults to 0.5.
#'
#' @return `select_coglasso()` returns an object of `S3` class `select_coglasso`
#'   containing the results of the
#'   selection procedure, built upon an object of `S3` class `coglasso`. Some
#'   output elements depend on the chosen model selection method. \cr 
#'   These elements are returned by all methods:
#' \itemize{
#'   \item ... are the same elements returned by [coglasso()].
#'   \item `sel_index_c`, `sel_index_lw` and `sel_index_lb` are the indexes of the
#'     final selected parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b}
#'     leading to the most stable sparse network.
#'   \item `sel_c`, `sel_lambda_w` and `sel_lambda_b` are the final selected
#'     parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} leading to the most
#'     stable sparse network.
#'   \item `sel_adj` is the adjacency matrix of the final selected network.
#'   \item `sel_density` is the density of the final selected network.
#'   \item `sel_icov` is the inverse covariance matrix of the final selected network.
#'   \item `sel_cov` optional, given only when `coglasso()` was called with 
#'   `cov_output = TRUE`. It is the covariance matrix associated with the final 
#'   selected network.
#'   \item `call` is the matched call.
#'   \item `method` is the chosen model selection method.
#' }
#'   These are the additional elements returned when choosing "xestars" or "xstars":
#' \itemize{
#'   \item `merge` is the "merged" adjacency matrix, the average of all the adjacency 
#'     matrices estimated across all the different subsamples for the selected 
#'     combination of \eqn{\lambda_w}, \eqn{\lambda_b}, and \eqn{c} values in the
#'     last path explored before convergence. Each entry is a measure of how 
#'     recurrent the corresponding edge is across the subsamples.
#'   \item `variability_lw`, `variability_lb` and `variability_c` are numeric vectors
#'     of as many items as the number of \eqn{\lambda_w}, \eqn{\lambda_b}, and 
#'     \eqn{c} values explored. Each item is the variability of the network 
#'     estimated for the corresponding hyperparameter value, keeping the other two 
#'     hyperparameters fixed to their selected value.
#'   \item `sel_variability` is the variability of the final selected network.
#' }
#'   These are the additional elements returned when choosing "ebic":
#' \itemize{
#'   \item `ebic_scores` is a numerical vector containing the eBIC scores for all the
#'     hyperparameter combination.
#' }
#' @export
#'
#' @examples
#' cg <- coglasso(multi_omics_sd_micro, p = c(4, 2), nlambda_w = 3, 
#'                nlambda_b = 3, nc = 3, verbose = FALSE)
#' # Using eXtended Efficient StARS, takes less than five seconds
#' sel_cg_xestars <- select_coglasso(cg, method = "xestars", verbose = FALSE)
#' \donttest{
#' # Using eXtended StARS, takes around a minute
#' sel_cg_xstars <- select_coglasso(cg, method = "xstars", verbose = FALSE)
#' }
#' # Using eBIC
#' sel_cg_ebic <- select_coglasso(cg, method = "ebic", verbose = FALSE)
#' 
select_coglasso <- function(coglasso_obj, method = "xestars", stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, max_iter = 10, old_sampling = FALSE, ebic_gamma = 0.5, verbose = TRUE) {
  call <- match.call()
  if ((method != "xestars") & (method != "xstars") & (method != "ebic")) {
    warning("Only available selection methods are \"xstars\", \"xestars\" and \"ebic\". Reverting to default selection method \"xestars\"")
    method <- "xestars"
  }

  if (method == "xstars") {
    coglasso_obj <- xstars(coglasso_obj = coglasso_obj, stars_thresh = stars_thresh, stars_subsample_ratio = stars_subsample_ratio, rep_num = rep_num, max_iter = max_iter, verbose = verbose)
  } else if (method == "xestars") {
    coglasso_obj <- xestars(coglasso_obj = coglasso_obj, stars_thresh = stars_thresh, stars_subsample_ratio = stars_subsample_ratio, rep_num = rep_num, max_iter = max_iter, old_sampling = old_sampling, verbose = verbose)
  } else if (method == "ebic") {
    if (verbose) {
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values with \"ebic\"....in progress"
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }
    
    n <- nrow(coglasso_obj$data)
    p_tot <- ncol(coglasso_obj$data)
    coglasso_obj$ebic_scores <- -n * coglasso_obj$loglik + log(n) * coglasso_obj$df + 4 * ebic_gamma * log(p_tot) * coglasso_obj$df
    
    if (verbose) {
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values with \"ebic\"....done"
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }

    sel_combination <- which.min(coglasso_obj$ebic_scores)
    sel_alpha <- coglasso_obj$hpars[sel_combination, 1]
    coglasso_obj$sel_index_c <- which(coglasso_obj$c == coglasso_obj$hpars[sel_combination, 4])
    coglasso_obj$sel_index_lw <- which(coglasso_obj$lambda_w == coglasso_obj$hpars[sel_combination, 2])
    coglasso_obj$sel_index_lb <- which(coglasso_obj$lambda_b == coglasso_obj$hpars[sel_combination, 3])
    coglasso_obj$sel_c <- coglasso_obj$c[coglasso_obj$sel_index_c]
    coglasso_obj$sel_lambda_w <- coglasso_obj$lambda_w[coglasso_obj$sel_index_lw]
    coglasso_obj$sel_lambda_b <- coglasso_obj$lambda_b[coglasso_obj$sel_index_lb]
    coglasso_obj$sel_adj <- coglasso_obj$path[[sel_combination]]
    coglasso_obj$sel_density <- coglasso_obj$density[sel_combination]
    coglasso_obj$sel_icov <- coglasso_obj$icov[[sel_combination]]
    if (!is.null(coglasso_obj$cov)) {
      coglasso_obj$sel_cov <- coglasso_obj$cov[[sel_combination]]
    }
  }
  coglasso_obj$call <- call
  coglasso_obj$method <- method
  
  class(coglasso_obj) <- "select_coglasso"
  
  return(coglasso_obj)
}


#' Print function for the S3 class `select_coglasso`
#'
#' Print information on the selected networks and the explored hyperparameters 
#' and see next suggested step
#'
#' @param x is the object of `S3` class `select_coglasso`.
#' @param ... system required, not used.
#' 
#' @noRd
#' @export
print.select_coglasso <- function(x, ...){
  sel_cg_name <- rlang::call_args(match.call())
  cat("Selected network estimated with collaborative graphical lasso\n\n")
  cat("The call was:\n", 
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("The model selection method was:\n", 
      x$method, "\n", sep = "")
  cat("The density of the selected network is:\n", 
      x$sel_density, "\n\n", sep = "")
  cat("Networks are made of ", x$D, " omics layers, for a total of ", ncol(x$data), " nodes\n", sep = "")
  mes <- paste(x$p[-length(x$p)], collapse = ", ")
  cat("For each layer they have: ", mes, " and ", x$p[length(x$p)], " nodes, respectively\n\n", sep = "")
  cat("The selected value for lambda within is:\n", 
      round(x$sel_lambda_w, 4), "\n", sep = "")
  cat("The selected value for lambda between is:\n", 
      round(x$sel_lambda_b, 4), "\n", sep = "")
  cat("The selected value for c is:\n", 
      round(x$sel_c, 4), "\n\n", sep = "")
  cat("The total number of hyperparameter combinations explored was:\n", 
      dim(x$hpars)[1], "\n", sep = "")
  cat("The values explored for lambda within were:\n", 
      paste(round(x$lambda_w, 4), collapse = ", "), "\n", sep = "")
  cat("The values explored for lambda between were:\n", 
      paste(round(x$lambda_b, 4), collapse = ", "), "\n", sep = "")
  cat("The values explored for c were:\n", paste(round(x$c, 4), collapse = ", "), "\n\n",
      sep = "")
  cat("Plot the selected network with:\nplot(", sel_cg_name[[1]], ")",
      sep = "")
}