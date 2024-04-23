#' Build multiple networks and select the best one from a multi-omics data set
#'
#' `bs()` wraps the two main functions of the package in a single one:
#' [coglasso()], to build multiple multi-omics networks, and [select_coglasso()]
#' to select the best one according to the chosen criterion. 
#' 
#' When using `bs()`, first, [coglasso()] estimates multiple multi-omics networks
#' with the algorithm *collaborative graphical lasso*, one for each combination
#' of input values for the hyperparameters \eqn{\lambda_w}, \eqn{\lambda_b} and 
#' \eqn{c}. Then, [select_coglasso()] selects the best combination of 
#' hyperparameters given to `coglasso()` according to the selected model 
#' selection method. The three availble options that can be set for the argument
#' `method` are  "xstars", "xestars" and "ebic". For more information on these 
#' selection methods, visit the help page of [select_coglasso()].
#'
#' @encoding UTF-8
#' @inherit coglasso
#' @inherit select_coglasso
#' @param verbose Print information regarding the network building and the 
#' network selection processes.
#' 
#' @return `bs()` returns an object of `S3` class `select_coglasso` containing 
#'   several elements. The most 
#'   important is arguably `sel_adj`, the adjacency matrix of the 
#'   selected network. Some output elements depend on the chosen model selection
#'   method.\cr
#'   These elements are always returned, and they are the result of network 
#'   estimation with [coglasso()]:
#' \itemize{
#'   \item `loglik` is a numerical vector containing the \eqn{log} likelihoods of all
#'   the estimated networks.
#'   \item `density` is a numerical vector containing a measure of the density of all
#'   the estimated networks.
#'   \item `df` is an integer vector containing the degrees of freedom of all the
#'   estimated networks.
#'   \item `convergence` is a binary vector containing whether a network was
#'   successfully estimated for the given combination of hyperparameters or not.
#'   \item `path` is a list containing the adjacency matrices of all the estimated
#'   networks.
#'   \item `icov` is a list containing the inverse covariance matrices of all the
#'   estimated networks.
#'   \item `nexploded` is the number of combinations of hyperparameters for which
#'   `coglasso()` failed to converge.
#'   \item `data` is the input multi-omics data set.
#'   \item `hpars` is the ordered table of all the combinations of hyperparameters
#'   given as input to `bs()`, with \eqn{\alpha(\lambda_w+\lambda_b)}
#'   being the key to sort rows.
#'   \item `lambda_w`, `lambda_b`, and `c` are numerical vectors with, 
#'   respectively, all the \eqn{\lambda_w}, \eqn{\lambda_b}, and \eqn{c} values
#'   `bs()` used.
#'   \item `pX` is the number of variables of the first data set.
#'   \item `cov` optional, returned when `cov_output` is TRUE, is a list containing
#'   the variance-covariance matrices of all the estimated networks.
#' }
#'   These elements are returned by all selection methods available:
#' \itemize{
#'   \item `sel_index_c`, `sel_index_lw` and `sel_index_lb` are the indexes of the
#'     final selected parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b}
#'     leading to the most stable sparse network.
#'   \item `sel_c`, `sel_lambda_w` and `sel_lambda_b` are the final selected
#'     parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} leading to the most
#'     stable sparse network.
#'   \item `sel_adj` is the adjacency matrix of the final selected network.
#'   \item `sel_density` is the density of the final selected network.
#'   \item `sel_icov` is the inverse covariance matrix of the final selected network.
#'   \item `call` is the matched call.
#'   \item `method` is the chosen model selection method.
#' }
#'   These are the additional elements returned when choosing "xestars":
#' \itemize{
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
#'   \item `merge_lw` and `merge_lb` are returned only if `light` is set to FALSE.
#'     They are lists with as many elements as the number of
#'     \eqn{c} parameters explored. Every element is a "merged" adjacency matrix,
#'     the average of all the adjacency matrices estimated  for those specific 
#'     \eqn{c} and the selected \eqn{\lambda_w} (or \eqn{\lambda_b}) values 
#'     across all the subsampling in the last path explored before convergence, 
#'     the one when the final combination of \eqn{\lambda_w} and \eqn{\lambda_b} 
#'     is selected for the given \eqn{c} value.
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
#' # Typical usage: give the input data set, set the value for `pX` and the 
#' # number of hyperparameters to explore (to choose how extensively to explore 
#' # the possible hyperparameters). Then, let the default behavior do the rest:
#' 
#' sel_mo_net <- bs(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3,
#'                  nc = 3, verbose = FALSE)
#' 
bs <- function(data, pX, lambda_w = NULL, lambda_b = NULL, c = NULL, 
               nlambda_w = NULL, nlambda_b = NULL, nc = NULL, 
               lambda_w_max = NULL, lambda_b_max = NULL, c_max = NULL, 
               lambda_w_min_ratio = NULL, lambda_b_min_ratio = NULL, 
               c_min_ratio = NULL, cov_output = FALSE, method = "xestars", 
               stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, 
               max_iter = 10, old_sampling = FALSE, light = TRUE, 
               ebic_gamma = 0.5, verbose = TRUE){
  
  call <- match.call()
  
  cg <- coglasso(data = data, pX = pX, lambda_w = lambda_w, lambda_b = lambda_b,
                 c = c, nlambda_w = nlambda_w, nlambda_b = nlambda_b, nc = nc, 
                 lambda_w_max = lambda_w_max, lambda_b_max = lambda_b_max, 
                 c_max = c_max, lambda_w_min_ratio = lambda_w_min_ratio, 
                 lambda_b_min_ratio = lambda_b_min_ratio, 
                 c_min_ratio = c_min_ratio, cov_output = cov_output, 
                 verbose = verbose)
  
  cg <- select_coglasso(cg, method = method, stars_thresh = stars_thresh, 
                        stars_subsample_ratio = stars_subsample_ratio, 
                        rep_num = rep_num, max_iter = max_iter, 
                        old_sampling = old_sampling, light = light, 
                        ebic_gamma = ebic_gamma, verbose = verbose)
  
  cg$call <- call
  
  cg
}