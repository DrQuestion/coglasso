#' Select the best `coglasso` network
#'
#' `select_coglasso()` selects the best combination of hyperparameters given to
#' `coglasso()` according to the selected model selection method. The two
#' availble options that can be set for the argument `method` are  "xstars" and
#' "ebic".
#'
#' `select_coglasso()` provides two model selection strategies: "xstars" and
#' "ebic":
#' \itemize{
#'   \item "xstars" uses *eXtended StARS* (*XStARS*) selecting the most stable, yet sparse
#'     network. Stability is computed upon network estimation from multiple subsamples of the
#'     multi-omics data set, allowing repetition. Subsamples are collected for a
#'     fixed amount of times (`rep_num`), and with a fixed proportion of the total
#'     number of samples (`stars_subsample_ratio`). See [stars_coglasso()] for more
#'     information on the methodology.\cr
#'   \item "ebic" uses the *extended Bayesian Information*
#'     *Criterion* (*eBIC*) selecting the network that minimizes it. `gamma` sets the
#'     wait given to the extended component, turning the model selection method to
#'     the standard *BIC* if set to 0.
#' }
#'
#' @inherit stars_coglasso
#' @param method The model selection method to select the best combination of
#'   hyperparameters. The available options are "xstars" and "eBIC". Defaults to
#'   "xstars".
#' @param ebic_gamma The \eqn{\gamma} tuning parameter for *eBIC* selection, to set
#'   between 0 and 1. When set to 0 one has the standard *BIC*. Defaults to 0.5.
#'
#' @return `select_coglasso()` returns a list containing the results of the
#'   selection procedure, built upon the list returned by `coglasso()`. Some
#'   output elements are there only when "xstars" is chosen as method. \cr These
#'   elements are returned by all methods:
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
#' }
#'   These are the additional elements returned when choosing "xstars":
#' \itemize{
#'   \item `merge_lw` and `merge_lb` are lists with as many elements as the number of
#'     \eqn{c} parameters explored. Every element is in turn a list of as many
#'     matrices as the number of \eqn{\lambda_w} (or \eqn{\lambda_b}) values
#'     explored. Each matrix is the "merged" adjacency matrix, the average of all
#'     the adjacency matrices estimated  for those specific \eqn{c} and
#'     \eqn{\lambda_w} (or \eqn{\lambda_b}) values across all the subsampling in
#'     the last path explored before convergence, the one when the final
#'     combination of \eqn{\lambda_w} and \eqn{\lambda_b} is selected for the
#'     given \eqn{c} value.
#'   \item `variability_lw` and `variability_lb` are lists with as many elements as
#'     the number of \eqn{c} parameters explored. Every element is a numeric
#'     vector of as many items as the number of \eqn{\lambda_w} (or
#'     \eqn{\lambda_b}) values explored. Each item is the variability of the
#'     network estimated for those specific \eqn{c} and \eqn{\lambda_w} (or
#'     \eqn{\lambda_b}) values in the last path explored before convergence, the
#'     one when the final combination of \eqn{\lambda_w} and \eqn{\lambda_b} is
#'     selected for the given \eqn{c} value.
#'   \item `opt_adj` is a list of the adjacency matrices finally selected for each
#'     \eqn{c} parameter explored.
#'   \item `opt_variability` is a numerical vector containing the variabilities
#'     associated to the adjacency matrices in `opt_adj`.
#'   \item `opt_index_lw` and `opt_index_lb` are integer vectors containing the
#'     index of the selected \eqn{\lambda_w}s (or \eqn{\lambda_b}s) for each
#'     \eqn{c} parameters explored.
#'   \item `opt_lambda_w` and `opt_lambda_b` are vectors containing the selected
#'     \eqn{\lambda_w}s (or \eqn{\lambda_b}s) for each \eqn{c} parameters
#'     explored.
#' }
#'   These are the additional elements returned when choosing "ebic":
#' \itemize{
#'   \item `ebic_scores` is a numerical vector containing the eBIC scores for all the
#'     hyperparameter combination.
#' }
#' @export
#'
#' @examples
#' cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3, verbose = FALSE)
#' \donttest{
#' # Using eXtended StARS, takes around 20 seconds
#' sel_cg_xstars <- select_coglasso(cg, method = "xstars", verbose = FALSE)
#' }
#' # Using eBIC
#' sel_cg_xstars <- select_coglasso(cg, method = "ebic", verbose = FALSE)
select_coglasso <- function(coglasso_obj, method = "xstars", stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, max_iter = 10, ebic_gamma = 0.5, verbose = TRUE) {
  if ((method != "xstars") & (method != "ebic")) {
    warning("Only available selection methods are \"xstars\" and \"ebic\". Reverting to default selection method \"xstars\"")
    method <- "xstars"
  }

  if (method == "xstars") {
    coglasso_obj <- stars_coglasso(coglasso_obj = coglasso_obj, stars_thresh = stars_thresh, stars_subsample_ratio = stars_subsample_ratio, rep_num = rep_num, max_iter = max_iter, verbose = verbose)
  } else if (method == "ebic") {
    if (verbose) {
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values with \"ebic\"....in progress"
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }
    
    n <- nrow(coglasso_obj$data)
    d <- ncol(coglasso_obj$data)
    coglasso_obj$ebic_scores <- -n * coglasso_obj$loglik + log(n) * coglasso_obj$df + 4 * ebic_gamma * log(d) * coglasso_obj$df
    
    if (verbose) {
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values with \"ebic\"....done"
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }

    sel_combination <- which.min(coglasso_obj$ebic_scores)
    sel_alpha <- coglasso_obj$hpars[sel_combination, 1]
    coglasso_obj$sel_index_c <- which(sort(unique(coglasso_obj$hpars[, 1])) == sel_alpha)
    coglasso_obj$sel_index_lw <- which(coglasso_obj$lambda_w == coglasso_obj$hpars[sel_combination, 2])
    coglasso_obj$sel_index_lb <- which(coglasso_obj$lambda_b == coglasso_obj$hpars[sel_combination, 3])
    coglasso_obj$sel_c <- coglasso_obj$c[coglasso_obj$sel_index_c]
    coglasso_obj$sel_lambda_w <- coglasso_obj$c[coglasso_obj$sel_index_lw]
    coglasso_obj$sel_lambda_b <- coglasso_obj$c[coglasso_obj$sel_index_lb]
    coglasso_obj$sel_adj <- coglasso_obj$path[[sel_combination]]
    coglasso_obj$sel_density <- coglasso_obj$density[sel_combination]
    coglasso_obj$sel_icov <- coglasso_obj$icov[[sel_combination]]
  }
  
  return(coglasso_obj)
}
