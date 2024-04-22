#' Estimate networks from a multi-omics data set
#'
#' `coglasso()` estimates multiple multi-omics networks with the algorithm
#' *collaborative graphical lasso*, one for each combination of input values for
#' the hyperparameters \eqn{\lambda_w}, \eqn{\lambda_b} and \eqn{c}.
#'
#' @encoding UTF-8
#' @param data The input multi-omics data set. Rows should be samples, columns
#'   should be variables. Variables should be grouped by their assay (e.g.
#'   transcripts first, then metabolites). `data` is a required parameter.
#' @param pX The number of variables of the first data set (e.g. the number of
#'   transcripts). `pX` is a required parameter.
#' @param lambda_w A vector of values for the parameter \eqn{\lambda_w}, the
#'   penalization parameter for the "within" interactions. Overrides
#'   `nlambda_w`.
#' @param lambda_b A vector of values for the parameter \eqn{\lambda_b}, the
#'   penalization parameter for the "between" interactions. Overrides
#'   `nlambda_b`.
#' @param c A vector of values for the parameter \eqn{c}, the weight given to
#'   collaboration. Overrides `nc`.
#' @param nlambda_w The number of requested \eqn{\lambda_w} parameters to
#'   explore. A sequence of size `nlambda_w` of \eqn{\lambda_w} parameters will
#'   be generated. Defaults to 8. Ignored when `lambda_w` is set by the user.
#' @param nlambda_b The number of requested \eqn{\lambda_b} parameters to
#'   explore. A sequence of size `nlambda_b` of \eqn{\lambda_b} parameters will
#'   be generated. Defaults to 8. Ignored when `lambda_b` is set by the user.
#' @param nc The number of requested \eqn{c} parameters to explore. A sequence
#'   of size `nc` of \eqn{c} parameters will be generated. Defaults to 8.
#'   Ignored when `c` is set by the user.
#' @param lambda_w_max The greatest generated \eqn{\lambda_w}. By default it is
#'   computed with a data-driven approach. Ignored when `lambda_w` is set by the
#'   user.
#' @param lambda_b_max The greatest generated \eqn{\lambda_b}. By default it is
#'   computed with a data-driven approach. Ignored when `lambda_b` is set by the
#'   user.
#' @param c_max The greatest generated \eqn{c}. Defaults to 10. Ignored when `c`
#'   is set by the user.
#' @param lambda_w_min_ratio The ratio of the smallest generated \eqn{\lambda_w}
#'   over the greatest generated \eqn{\lambda_w}. Defaults to 0.1. Ignored when
#'   `lambda_w` is set by the user.
#' @param lambda_b_min_ratio The ratio of the smallest generated \eqn{\lambda_b}
#'   over the greatest generated \eqn{\lambda_b}. Defaults to 0.1. Ignored when
#'   `lambda_b` is set by the user.
#' @param c_min_ratio The ratio of the smallest generated \eqn{c} over the
#'   greatest generated \eqn{c}. Defaults to 0.1. Ignored when `c` is set by the
#'   user.
#' @param cov_output Add the estimated variance-covariance matrix to the output.
#' @param verbose Print information regarding current `coglasso` run on the
#'   console.
#'
#' @return `coglasso()` returns a list containing several elements:
#' * `loglik` is a numerical vector containing the \eqn{log} likelihoods of all
#'   the estimated networks.
#' * `density` is a numerical vector containing a measure of the density of all
#'   the estimated networks.
#' * `df` is an integer vector containing the degrees of freedom of all the
#'   estimated networks.
#' * `convergence` is a binary vector containing whether a network was
#'   successfully estimated for the given combination of hyperparameters or not.
#' * `path` is a list containing the adjacency matrices of all the estimated
#'   networks.
#' * `icov` is a list containing the inverse covariance matrices of all the
#'   estimated networks.
#' * `nexploded` is the number of combinations of hyperparameters for which
#'   `coglasso()` failed to converge.
#' * `data` is the input multi-omics data set.
#' * `hpars` is the ordered table of all the combinations of hyperparameters
#'   given as input to `coglasso()`, with \eqn{\alpha(\lambda_w+\lambda_b)}
#'   being the key to sort rows.
#' * `lambda_w` is a numerical vector with all the \eqn{\lambda_w} values `coglasso()`
#'   used.
#' * `lambda_b` is a numerical vector with all the \eqn{\lambda_b} values `coglasso()`
#'   used.
#' * `c` is a numerical vector with all the \eqn{c} values `coglasso()`
#'   used.
#' * `pX` is the number of variables of the first data set.
#' * `cov` optional, returned when `cov_output` is TRUE, is a list containing
#'   the variance-covariance matrices of all the estimated networks.
#'
#' @export
#'
#' @examples
#' # Typical usage: set the number of hyperparameters to explore
#' cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3, verbose = FALSE)
#' \donttest{
#' # Model selection using eXtended Efficient StARS, takes less than five seconds
#' sel_cg_xestars <- select_coglasso(cg, method = "xestars", verbose = FALSE)
#' }
#' 
coglasso <- function(data, pX, lambda_w = NULL, lambda_b = NULL, c = NULL, nlambda_w = NULL, nlambda_b = NULL, nc = NULL, lambda_w_max = NULL, lambda_b_max = NULL, c_max = NULL, lambda_w_min_ratio = NULL, lambda_b_min_ratio = NULL, c_min_ratio = NULL, cov_output = FALSE, verbose = TRUE) {
  original_data <- data
  data <- as.matrix(data)
  gcinfo(FALSE)
  n <- nrow(data)
  p <- ncol(data)
  cov.input <- isSymmetric(data)
  if (cov.input) {
    if (verbose) cat("The input is identified as the covariance matrix.\n")
    S <- data
  } else {
    data <- scale(data)
    S <- cor(data)
  }
  rm(data)
  gc()

  hpars <- gen_hpars(
    S = S, p = p, lambda_w = lambda_w, lambda_b = lambda_b, c = c, nlambda_w = nlambda_w, nlambda_b = nlambda_b, nc = nc, lambda_w_max = lambda_w_max, lambda_b_max = lambda_b_max, c_max = c_max,
    lambda_w_min_ratio = lambda_w_min_ratio, lambda_b_min_ratio = lambda_b_min_ratio, c_min_ratio = c_min_ratio
  )

  cg <- co_glasso(S, pX, hpars[[1]], FALSE, verbose, cov_output)

  cg$data <- original_data
  cg$hpars <- hpars[[1]]
  cg$lambda_w <- hpars[[2]]
  cg$lambda_b <- hpars[[3]]
  cg$c <- hpars[[4]]
  cg$pX <- pX

  cg
}

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
#' @return `bs()` returns a list containing several elements. The most 
#'   important is arguably `sel_adj`, that is the adjacency matrix of the 
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
#' sel_mo_net <- bs(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3, verbose = FALSE)
#' 
bs <- function(data, pX, lambda_w = NULL, lambda_b = NULL, c = NULL, 
               nlambda_w = NULL, nlambda_b = NULL, nc = NULL, 
               lambda_w_max = NULL, lambda_b_max = NULL, c_max = NULL, 
               lambda_w_min_ratio = NULL, lambda_b_min_ratio = NULL, 
               c_min_ratio = NULL, cov_output = FALSE, method = "xestars", 
               stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, 
               max_iter = 10, old_sampling = FALSE, light = TRUE, 
               ebic_gamma = 0.5, verbose = TRUE){

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
  
  cg
}

#' Generate combinations of hyperparameters for `coglasso()`
#'
#' `gen_hpars()` generates an ordered table of all the combinations of the hyperparameters given as input to `coglasso()`
#'
#' @inherit coglasso
#' @param S The empirical Pearson's correlation matrix of the multi-omics data set.
#' @param p The total number of variables in the multi-omics data set.
#' @param getc Parameter for testing purposes. If TRUE, the returned table contains a column with *c* values rather than *\alpha* values.
#'
#' @return
#' `gen_hpars()` returns a list of four elements:
#' * A table of all the combinations of hyperparameters given as input to `coglasso()`, with \eqn{\alpha(\lambda_w+\lambda_b)} being the key to sort rows.
#' * A numerical vector with all the generated \eqn{\lambda_w}.
#' * A numerical vector with all the generated \eqn{\lambda_b}.
#' * A numerical vector with all the generated \eqn{c}.
#' 
#' @noRd
gen_hpars <- function(S = NULL, p = NULL, lambda_w = NULL, lambda_b = NULL, c = NULL,
                      nlambda_w = NULL, nlambda_b = NULL, nc = NULL, 
                      lambda_w_max = NULL, lambda_b_max = NULL, c_max = NULL, 
                      lambda_w_min_ratio = NULL, lambda_b_min_ratio = NULL,
                      c_min_ratio = NULL, getc = NULL) {
  if (!is.null(lambda_w)) nlambda_w <- length(lambda_w)
  if (!is.null(lambda_b)) nlambda_b <- length(lambda_b)
  if (!is.null(c)) nalpha <- length(c)

  if (is.null(lambda_w)) {
    if (is.null(nlambda_w)) {
      nlambda_w <- 8
    }
    if (is.null(lambda_w_min_ratio)) {
      lambda_w_min_ratio <- 0.1
    }
    if (is.null(lambda_w_max)) {
      lambda_w_max <- max(max(S - diag(p)), -min(S - diag(p)))
    }
    lambda_w_min <- lambda_w_min_ratio * lambda_w_max
    lambda_w <- exp(seq(log(lambda_w_max), log(lambda_w_min), length = nlambda_w))
  }
  if (is.null(lambda_b)) {
    if (is.null(nlambda_b)) {
      nlambda_b <- 8
    }
    if (is.null(lambda_b_min_ratio)) {
      lambda_b_min_ratio <- 0.1
    }
    if (is.null(lambda_b_max)) {
      lambda_b_max <- max(max(S - diag(p)), -min(S - diag(p)))
    }
    lambda_b_min <- lambda_b_min_ratio * lambda_b_max
    lambda_b <- exp(seq(log(lambda_b_max), log(lambda_b_min), length = nlambda_b))
  }
  if (is.null(c)) {
    if (is.null(nc)) {
      nc <- 8
    }
    if (is.null(c_min_ratio)) {
      c_min_ratio <- 0.1
    }
    if (is.null(c_max)) {
      c_max <- 10
    }
    c.min <- c_min_ratio * c_max
    c <- exp(seq(log(c_max), log(c.min), length = nc))
  }

  hpars <- vector(mode = "list", length = 4)

  if (is.null(getc)) {
    alpha <- 1 / (1 + c)
    hpars_table <- expand.grid(alpha, lambda_w, lambda_b)
    key <- hpars_table[, 1] * (hpars_table[, 2] + hpars_table[, 3])
    hpars_table <- cbind(hpars_table, key)
    hpars_table <- hpars_table[order(hpars_table$key, decreasing = T), ]
  } else {
    hpars_table <- expand.grid(c, lambda_w, lambda_b)
    key <- hpars_table[, 1] * (hpars_table[, 2] + hpars_table[, 3])
    hpars_table <- cbind(hpars_table, key)
    hpars_table <- hpars_table[order(hpars_table$key, decreasing = F), ]
  }

  hpars_table <- hpars_table[, c(1, 2, 3)]
  hpars[[1]] <- as.matrix(hpars_table)
  if (is.null(getc)) {
    colnames(hpars[[1]]) <- c("alpha", "lambda_w", "lambda_b")
  } else {
    colnames(hpars[[1]]) <- c("c", "lambda_w", "lambda_b")
  }
  rownames(hpars[[1]]) <- seq(1:nrow(hpars[[1]]))
  hpars[[2]] <- lambda_w
  hpars[[3]] <- lambda_b
  hpars[[4]] <- c

  hpars
}
