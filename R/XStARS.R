#' Stability selection of the best `coglasso` network
#' 
#' @description 
#' `xstars()` selects the combination of hyperparameters given to
#' `coglasso()` yielding the most stable, yet sparse network. Stability is
#' computed upon network estimation from multiple subsamples of the multi-omics data set,
#' allowing repetition. Subsamples are collected for a fixed amount of times
#' (`rep_num`), and with a fixed proportion of the total number of samples
#' (`stars_subsample_ratio`).
#'
#' @details 
#' *eXtended StARS* (*XStARS*) is an adaptation for *collaborative graphical regression* of the method
#' published by Liu, H. *et al.* (2010): Stability Approach to Regularization
#' Selection (StARS). *StARS* was developed for network estimation regulated by
#' a single penalty parameter, while *collaborative graphical lasso* needs to
#' explore three different hyperparameters. In particular, two of these are
#' penalty parameters with a direct influence on network sparsity, hence on
#' stability. For every \eqn{c} parameter, `xstars()` explores one of
#' the two penalty parameters (\eqn{\lambda_w} or \eqn{\lambda_b}), keeping the other one
#' fixed at its previous best estimate, using the normal, one-dimentional
#' *StARS* approach, until finding the best couple. It then selects the \eqn{c}
#' parameter for which the best (\eqn{\lambda_w}, \eqn{\lambda_b}) couple yielded the most
#' stable, yet sparse network.
#'
#' @param coglasso_obj The object of `S3` class `coglasso` returned by `coglasso()`.
#' @param stars_thresh The threshold set for variability of the explored
#'   networks at each iteration of the algorithm. The \eqn{\lambda_w} or the \eqn{\lambda_b}
#'   associated to the most stable network before the threshold is overcome is
#'   selected.
#' @param stars_subsample_ratio The proportion of samples in the multi-omics
#'   data set to be randomly subsampled to estimate the variability of the
#'   network under the given hyperparameters setting. Defaults to 80% when the
#'   number of samples is smaller than 144, otherwise it defaults to
#'   \eqn{\frac{10}{n}\sqrt{n}}.
#' @param rep_num The amount of subsamples of the multi-omics data set used to
#'   estimate the variability of the network under the given hyperparameters
#'   setting. Defaults to 20.
#' @param max_iter The greatest number of times the algorithm is allowed to
#'   choose a new best \eqn{\lambda_w}. Defaults to 10.
#' @param verbose Print information regarding the progress of the selection
#'   procedure on the console.
#'
#' @return `xstars()` returns an object of `S3` class `select_coglasso`
#'   containing the results of the
#'   selection procedure, built upon the object of `S3` class `coglasso` returned by `coglasso()`.
#' * ... are the same elements returned by [coglasso()].
#' * `merge_lw` and `merge_lb` are lists with as many elements as the number of
#'   \eqn{c} parameters explored. Every element is in turn a list of as many
#'   matrices as the number of \eqn{\lambda_w} (or \eqn{\lambda_b}) values explored. Each
#'   matrix is the "merged" adjacency matrix, the average of all the adjacency
#'   matrices estimated  for those specific \eqn{c} and \eqn{\lambda_w} (or \eqn{\lambda_b})
#'   values across all the subsampling in the last path explored before
#'   convergence, the one when the final combination of \eqn{\lambda_w} and \eqn{\lambda_b}
#'   is selected for the given \eqn{c} value.
#' * `variability_lw` and `variability_lb` are lists with as many elements as
#'  the number of \eqn{c} parameters explored. Every element is a numeric vector
#'  of as many items as the number of \eqn{\lambda_w} (or \eqn{\lambda_b}) values explored.
#'  Each item is the variability of the network estimated for those specific
#'  \eqn{c} and \eqn{\lambda_w} (or \eqn{\lambda_b}) values in the last path explored before
#'  convergence, the one when the final combination of \eqn{\lambda_w} and \eqn{\lambda_b}
#'  is selected for the given \eqn{c} value.
#' * `opt_adj` is a list of the adjacency matrices finally selected for each
#'  \eqn{c} parameter explored.
#' * `opt_variability` is a numerical vector containing the variabilities
#'  associated to the adjacency matrices in `opt_adj`.
#' * `opt_index_lw` and `opt_index_lb` are integer vectors containing the
#'  index of the selected \eqn{\lambda_w}s (or \eqn{\lambda_b}s) for each \eqn{c} parameters
#'  explored.
#' * `opt_lambda_w` and `opt_lambda_b` are vectors containing the selected
#'  \eqn{\lambda_w}s (or \eqn{\lambda_b}s) for each \eqn{c} parameters explored.
#' * `sel_index_c`, `sel_index_lw` and `sel_index_lb` are the indexes of the
#'  final selected parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} leading to the
#'  most stable sparse network.
#' * `sel_c`, `sel_lambda_w` and `sel_lambda_b` are the final selected
#'  parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} leading to the most stable
#'  sparse network.
#' * `sel_adj` is the adjacency matrix of the final selected network.
#' * `sel_density` is the density of the final selected network.
#' * `sel_icov` is the inverse covariance matrix of the final selected network.
#' * `call` is the matched call.
#' * `method` is the chosen model selection method. Here, it is "xstars".
#'
#' @export
#'
#' @examples
#' cg <- coglasso(multi_omics_sd_micro, p = c(4, 2), nlambda_w = 3, 
#'                nlambda_b = 3, nc = 3, verbose = FALSE)
#' \donttest{
#' # Takes around one minute
#' sel_cg <- xstars(cg, verbose = FALSE)
#' }
xstars <- function(coglasso_obj, stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, max_iter = 10, verbose = TRUE) {
  call <- match.call()
  
  n <- nrow(coglasso_obj$data)
  p_tot <- ncol(coglasso_obj$data)
  p <- coglasso_obj$p
  D <- coglasso_obj$D
  n_lambda_w <- length(coglasso_obj$lambda_w)
  n_lambda_b <- length(coglasso_obj$lambda_b)
  n_c <- length(coglasso_obj$c)
  
  if (is.null(stars_subsample_ratio)) {
    if (n > 144) stars_subsample_ratio <- 10 * sqrt(n) / n
    if (n <= 144) stars_subsample_ratio <- 0.8
  }
  
  coglasso_obj$merge_lw <- vector(mode = "list", length = n_c)
  coglasso_obj$merge_lb <- vector(mode = "list", length = n_c)
  coglasso_obj$variability_lw <- vector(mode = "list", length = n_c)
  coglasso_obj$variability_lb <- vector(mode = "list", length = n_c)
  coglasso_obj$opt_adj <- vector(mode = "list", length = n_c)
  coglasso_obj$opt_variability <- rep(0, n_c)
  coglasso_obj$opt_index_lw <- rep(0, n_c)
  coglasso_obj$opt_index_lb <- rep(0, n_c)
  coglasso_obj$opt_lambda_w <- rep(0, n_c)
  coglasso_obj$opt_lambda_b <- rep(0, n_c)
  
  for (i in 1:n_c) {
    if (verbose) {
      mes <- paste(c("Selecting best Lambda_w/Lambda_b combination for all c values with \"xstars\"....in progress:", floor(100 * i / n_c), "%"), collapse = "")
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }
    
    c <- coglasso_obj$c[i]
    alpha <- 1 / (c * (D - 1) + 1)
    
    lw_sel <- -1
    lb_sel <- coglasso_obj$lambda_b[n_lambda_b]
    num_iter <- 0
    converged <- FALSE
    select_lw <- TRUE
    
    while (!converged & num_iter < max_iter) {
      if (select_lw) {
        coglasso_obj$merge_lw[[i]] <- list()
        for (j in 1:n_lambda_w) coglasso_obj$merge_lw[[i]][[j]] <- Matrix::Matrix(0, p_tot, p_tot)
        
        real_rep.num <- rep(0, n_lambda_w)
        
        for (j in 1:rep_num)
        {
          if (verbose) {
            mes <- paste(c("Conducting Subsampling Lambda_w....in progress:", floor(100 * j / rep_num), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
          }
          ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
          
          corr_matrix <- cor(scale(coglasso_obj$data[ind.sample, ]))
          hpars <- matrix(c(rep(alpha, n_lambda_w), coglasso_obj$lambda_w, rep(lb_sel, n_lambda_w), rep(c, n_lambda_w)), nrow = n_lambda_w, ncol = 4)
          tmp <- co_glasso_D(corr_matrix, p, hpars, FALSE, FALSE, FALSE)
          
          convergence <- tmp$convergence
          tmp <- tmp$path
          
          for (k in 1:n_lambda_w) {
            if (convergence[k] == 1) {
              real_rep.num[k] <- real_rep.num[k] + 1
              coglasso_obj$merge_lw[[i]][[k]] <- coglasso_obj$merge_lw[[i]][[k]] + tmp[[k]]
            }
          }
          rm(ind.sample, tmp)
          gc()
        }
        
        coglasso_obj$variability_lw[[i]] <- rep(0, n_lambda_w)
        for (j in 1:n_lambda_w) {
          coglasso_obj$merge_lw[[i]][[j]] <- coglasso_obj$merge_lw[[i]][[j]] / real_rep.num[j]
          coglasso_obj$variability_lw[[i]][j] <- 4 * sum(coglasso_obj$merge_lw[[i]][[j]] *
                                                           (1 - coglasso_obj$merge_lw[[i]][[j]])) / (p_tot * (p_tot - 1))
        }
        
        index_lw <- max(which.max(coglasso_obj$variability_lw[[i]] >=
                                    stars_thresh)[1] - 1, 1)
        coglasso_obj$opt_variability[i] <- coglasso_obj$variability_lw[[i]][index_lw]
        tmp_lw <- coglasso_obj$lambda_w[[index_lw]]
        if (tmp_lw == lw_sel) converged <- TRUE
        
        lw_sel <- tmp_lw
        
        if (verbose) {
          mes <- "Conducting Subsampling Lambda_w....done.                 "
          cat(mes, "\r")
          cat("\n")
          if (converged) {
            mes <- "Reached convergence.                 "
            cat(mes, "\r")
            cat("\n")
          }
          flush.console()
        }
      } else {
        coglasso_obj$merge_lb[[i]] <- list()
        for (j in 1:n_lambda_b) coglasso_obj$merge_lb[[i]][[j]] <- Matrix::Matrix(0, p_tot, p_tot)
        
        real_rep.num <- rep(0, n_lambda_b)
        
        for (j in 1:rep_num)
        {
          if (verbose) {
            mes <- paste(c("Conducting Subsampling Lambda_b....in progress:", floor(100 * j / rep_num), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
          }
          ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
          
          corr_matrix <- cor(scale(coglasso_obj$data[ind.sample, ]))
          hpars <- matrix(c(rep(alpha, n_lambda_b), coglasso_obj$lambda_b, rep(lw_sel, n_lambda_b), rep(c, n_lambda_b)), nrow = n_lambda_b, ncol = 4)
          tmp <- co_glasso_D(corr_matrix, p, hpars, FALSE, FALSE, FALSE)
          
          convergence <- tmp$convergence
          tmp <- tmp$path
          
          for (k in 1:n_lambda_b) {
            if (convergence[k] == 1) {
              real_rep.num[k] <- real_rep.num[k] + 1
              coglasso_obj$merge_lb[[i]][[k]] <- coglasso_obj$merge_lb[[i]][[k]] + tmp[[k]]
            }
          }
          rm(ind.sample, tmp)
          gc()
        }
        
        coglasso_obj$variability_lb[[i]] <- rep(0, n_lambda_b)
        for (j in 1:n_lambda_b) {
          coglasso_obj$merge_lb[[i]][[j]] <- coglasso_obj$merge_lb[[i]][[j]] / real_rep.num[j]
          coglasso_obj$variability_lb[[i]][j] <- 4 * sum(coglasso_obj$merge_lb[[i]][[j]] *
                                                           (1 - coglasso_obj$merge_lb[[i]][[j]])) / (p_tot * (p_tot - 1))
        }
        
        index_lb <- max(which.max(coglasso_obj$variability_lb[[i]] >=
                                    stars_thresh)[1] - 1, 1)
        coglasso_obj$opt_variability[i] <- coglasso_obj$variability_lb[[i]][index_lb]
        tmp_lb <- coglasso_obj$lambda_b[[index_lb]]
        if (tmp_lb == lb_sel) converged <- TRUE
        
        lb_sel <- tmp_lb
        
        if (verbose) {
          mes <- "Conducting Subsampling Lambda_b....done.                 "
          cat(mes, "\r")
          cat("\n")
          if (converged) {
            mes <- "Reached convergence.                 "
            cat(mes, "\r")
            cat("\n")
          }
          flush.console()
        }
      }
      select_lw <- !select_lw
      
      num_iter <- num_iter + 0.5
    }
    if (num_iter == max_iter) {
      if (verbose) {
        mes <- "Reached max Iterations.                 "
        cat(mes, "\r")
        cat("\n")
        flush.console()
      }
    }
    
    
    if (verbose) {
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values with \"xstars\"....done"
      cat(mes, "\r")
      flush.console()
    }
    
    coglasso_obj$opt_index_lw[i] <- index_lw
    coglasso_obj$opt_index_lb[i] <- index_lb
    coglasso_obj$opt_lambda_w[i] <- coglasso_obj$lambda_w[index_lw]
    coglasso_obj$opt_lambda_b[i] <- coglasso_obj$lambda_b[index_lb]
    opt_hpars_combination <- which(coglasso_obj$hpars[, 1] == alpha &
                                     coglasso_obj$hpars[, 2] == coglasso_obj$opt_lambda_w[i] &
                                     coglasso_obj$hpars[, 3] == coglasso_obj$opt_lambda_b[i] &
                                     coglasso_obj$hpars[, 4] == c)
    coglasso_obj$opt_adj[[i]] <- coglasso_obj$path[[opt_hpars_combination]]
    colnames(coglasso_obj$opt_adj[[i]]) <- colnames(coglasso_obj$data)
    row.names(coglasso_obj$opt_adj[[i]]) <- colnames(coglasso_obj$data)
  }
  
  if (all(is.na(coglasso_obj$opt_variability))) {
    warning("coglasso did not converge for any c parameter. Will select the highest c.")
  }
  
  coglasso_obj$sel_index_c <- min(which.min(coglasso_obj$opt_variability), n_c)
  coglasso_obj$sel_index_lw <- coglasso_obj$opt_index_lw[coglasso_obj$sel_index_c]
  coglasso_obj$sel_index_lb <- coglasso_obj$opt_index_lb[coglasso_obj$sel_index_c]
  coglasso_obj$sel_lambda_w <- coglasso_obj$lambda_w[coglasso_obj$sel_index_lw]
  coglasso_obj$sel_lambda_b <- coglasso_obj$lambda_b[coglasso_obj$sel_index_lb]
  coglasso_obj$sel_c <- coglasso_obj$c[coglasso_obj$sel_index_c]
  coglasso_obj$sel_adj <- coglasso_obj$opt_adj[[coglasso_obj$sel_index_c]]
  sel_hpars_combination <- which(coglasso_obj$hpars[, 1] == 1 / (coglasso_obj$sel_c * (D - 1) +1) &
                                   coglasso_obj$hpars[, 2] == coglasso_obj$sel_lambda_w &
                                   coglasso_obj$hpars[, 3] == coglasso_obj$sel_lambda_b &
                                   coglasso_obj$hpars[, 4] == coglasso_obj$sel_c)
  coglasso_obj$sel_density <- coglasso_obj$density[sel_hpars_combination]
  coglasso_obj$sel_icov <- coglasso_obj$icov[[sel_hpars_combination]]
  if (!is.null(coglasso_obj$cov)) {
    coglasso_obj$sel_cov <- coglasso_obj$cov[[sel_hpars_combination]]
  }
  coglasso_obj$method <- "xstars"
  coglasso_obj$call <- call
  
  class(coglasso_obj) <- "select_coglasso"
  
  return(coglasso_obj)
}

#' Efficient stability selection of the best `coglasso` network
#'
#' @description `xestars()` provides a more efficient and lighter implementation
#' than `xstars()` to select the combination of hyperparameters given to
#' `coglasso()` yielding the most stable, yet sparse network. Stability is
#' computed upon network estimation from multiple subsamples of the multi-omics
#' data set, allowing repetition. Subsamples are collected for a fixed amount of
#' times (`rep_num`), and with a fixed proportion of the total number of samples
#' (`stars_subsample_ratio`).
#'
#' @details
#' *eXtended Efficient StARS* (*XEStARS*) is a more efficient and memory-light version of
#' *XStARS*, the adaptation for *collaborative graphical regression* of the method
#' published by Liu, H. *et al.* (2010): Stability Approach to Regularization
#' Selection (StARS). *StARS* was developed for network estimation regulated by
#' a single penalty parameter, while *collaborative graphical lasso* needs to
#' explore three different hyperparameters. In particular, two of these are
#' penalty parameters with a direct influence on network sparsity, hence on
#' stability. For every \eqn{c} parameter, `xestars()` explores one of the two
#' penalty parameters (\eqn{\lambda_w} or \eqn{\lambda_b}), keeping the other
#' one fixed at its previous best estimate, using the normal, one-dimentional
#' *StARS* approach, until finding the best couple. What makes it more efficient
#' than `xstars()` is that the stability check that in the original algorithm 
#' (even in the original *StARS*) is performed for every \eqn{\lambda_w} or 
#' \eqn{\lambda_b} value, is implemented here as a *stopping criterion*. This 
#' reduces sensibly the number of iterations before convergence. It then selects
#' the \eqn{c} parameter for which the best (\eqn{\lambda_w}, \eqn{\lambda_b})
#' couple yielded the most stable, yet sparse network. \cr
#' The original *XStARS* computes a new subsampling for every time the algorithm
#' switches from optimizing the two\eqn{\lambda_w} and \eqn{\lambda_b}, and for 
#' every \eqn{c}. This does not allow to compare the hyperparameters on an equal
#' ground, and can slow the selection down with bigger data set or a larger
#' hyperparameter space. To allow a fairer (and faster) comparison among 
#' different optimizations, the `old_sampling` parameter has been implemented. 
#' If set to TRUE, the subsampling is the same one `xstars()` would perform. 
#' Otherwise the subsampling is performed at the beginning of the algorithm once
#' and for all its iterations. \cr
#' To allow `xestars()` to be more memory light, the `light` parameter has been 
#' implemented. If set to TRUE and the "merged" matrixes traditionally returned 
#' by both *StARS* and *XStARS* are not returned.
#' 
#' @param coglasso_obj The object of `S3` class `coglasso` returned by `coglasso()`.
#' @param stars_thresh The threshold set for variability of the explored
#'   networks at each iteration of the algorithm. The \eqn{\lambda_w} or the
#'   \eqn{\lambda_b} associated to the most stable network before the threshold
#'   is overcome is selected.
#' @param stars_subsample_ratio The proportion of samples in the multi-omics
#'   data set to be randomly subsampled to estimate the variability of the
#'   network under the given hyperparameters setting. Defaults to 80% when the
#'   number of samples is smaller than 144, otherwise it defaults to
#'   \eqn{\frac{10}{n}\sqrt{n}}.
#' @param rep_num The amount of subsamples of the multi-omics data set used to
#'   estimate the variability of the network under the given hyperparameters
#'   setting. Defaults to 20.
#' @param max_iter The greatest number of times the algorithm is allowed to
#'   choose a new best \eqn{\lambda_w}. Defaults to 10.
#' @param old_sampling Perform the same subsampling `xstars()` would if set to 
#'   TRUE. Makes a difference with bigger data sets, where computing
#'   a correlation matrix could take significantly longer. Defaults to FALSE.
#' @param light Do not store the "merged" matrixes recording average variability
#'   of each edge, making the algorithm more memory efficient, if set to TRUE. 
#'   Defaults to TRUE.
#' @param verbose Print information regarding the progress of the selection
#'   procedure on the console.
#'
#' @return `xestars()` returns an object of `S3` class `select_coglasso` 
#'   containing the results of the selection
#'   procedure, built upon the object of `S3` class `coglasso` returned by `coglasso()`.
#' * ... are the same elements returned by [coglasso()].
#' * `opt_adj` is a list of the adjacency matrices finally selected for each
#'   \eqn{c} parameter explored.
#' * `opt_variability` is a numerical vector containing the variabilities
#'   associated to the adjacency matrices in `opt_adj`.
#' * `opt_index_lw` and `opt_index_lb` are integer vectors containing the
#'   index of the selected \eqn{\lambda_w}s (or \eqn{\lambda_b}s) for each
#'   \eqn{c} parameters explored.
#' * `opt_lambda_w` and `opt_lambda_b` are vectors containing the selected
#'   \eqn{\lambda_w}s (or \eqn{\lambda_b}s) for each \eqn{c} parameters
#'   explored.
#' * `sel_index_c`, `sel_index_lw` and `sel_index_lb` are the indexes of the
#'   final selected parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b}
#'   leading to the most stable sparse network.
#' * `sel_c`, `sel_lambda_w` and `sel_lambda_b` are the final selected
#'   parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} leading to the most
#'   stable sparse network.
#' * `sel_adj` is the adjacency matrix of the final selected network.
#' * `sel_density` is the density of the final selected network.
#' * `sel_icov` is the inverse covariance matrix of the final selected network.
#' * `call` is the matched call.
#' * `method` is the chosen model selection method. Here, it is "xestars".
#' * `merge_lw` and `merge_lb` are returned only if `light` is set to FALSE.
#'   They are lists with as many elements as the number of
#'   \eqn{c} parameters explored. Every element is a "merged" adjacency matrix,
#'   the average of all the adjacency matrices estimated  for those specific 
#'   \eqn{c} and the selected \eqn{\lambda_w} (or \eqn{\lambda_b}) values 
#'   across all the subsampling in the last path explored before convergence, 
#'   the one when the final combination of \eqn{\lambda_w} and \eqn{\lambda_b} 
#'   is selected for the given \eqn{c} value.
#'
#' @export
#'
#' @examples
#' cg <- coglasso(multi_omics_sd_micro, p = c(4, 2), nlambda_w = 3, 
#'                nlambda_b = 3, nc = 3, verbose = FALSE)
#' \donttest{
#' # Takes less than five seconds
#' sel_cg <- xestars(cg, verbose = FALSE)
#' }
xestars <- function(coglasso_obj, stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, max_iter = 10, old_sampling = FALSE, light = TRUE, verbose = TRUE) {
  call <- match.call()
  
  n <- nrow(coglasso_obj$data)
  p_tot <- ncol(coglasso_obj$data)
  p <- coglasso_obj$p
  D <- coglasso_obj$D
  n_lambda_w <- length(coglasso_obj$lambda_w)
  n_lambda_b <- length(coglasso_obj$lambda_b)
  n_c <- length(coglasso_obj$c)
  
  if (is.null(stars_subsample_ratio)) {
    if (n > 144) stars_subsample_ratio <- 10 * sqrt(n) / n
    if (n <= 144) stars_subsample_ratio <- 0.8
  }
  
  if (!light) {
    coglasso_obj$merge_lw <- vector(mode = "list", length = n_c)
    coglasso_obj$merge_lb <- vector(mode = "list", length = n_c)
  }
  
  coglasso_obj$opt_adj <- vector(mode = "list", length = n_c)
  coglasso_obj$opt_variability <- rep(0, n_c)
  coglasso_obj$opt_index_lw <- rep(0, n_c)
  coglasso_obj$opt_index_lb <- rep(0, n_c)
  coglasso_obj$opt_lambda_w <- rep(0, n_c)
  coglasso_obj$opt_lambda_b <- rep(0, n_c)
  
  # Create the same samples of indexes for the whole set of iterations,
  # cancelled when `old_sampling` is true
  if (!old_sampling) {
    corr_matrixes <- vector(mode = "list", length = rep_num)
    for (j in 1:rep_num) {
      ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
      corr_matrixes[[j]] <- cor(scale(coglasso_obj$data[ind.sample, ]))
    }
  }
  
  for (i in 1:n_c) {
    if (verbose) {
      mes <- paste(c("Selecting best Lambda_w/Lambda_b combination for all c values with \"xestars\"....in progress:", floor(100 * i / n_c), "%"), collapse = "")
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }
    
    c <- coglasso_obj$c[i]
    alpha <- 1 / (c * (D - 1) + 1)
    
    lw_sel <- -1
    lb_sel <- coglasso_obj$lambda_b[n_lambda_b]
    num_iter <- 0
    converged <- FALSE
    select_lw <- TRUE
    
    while (!converged & num_iter < max_iter) {
      if (select_lw) {
        # Consider breaking in a one-dimensional sub-function returning a lambda estars
        
        if (old_sampling) {
          corr_matrixes <- vector(mode = "list", length = rep_num)
          for (j in 1:rep_num) {
            ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
            corr_matrixes[[j]] <- cor(scale(coglasso_obj$data[ind.sample, ]))
          }
        }
        
        for (j in 1:n_lambda_w) {
          if (verbose) {
            mes <- paste(c("Conducting Subsampling Lambda_w....in progress:", floor(100 * j / n_lambda_w), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
          }
          real_rep.num <- 0
          merge_tmp <- rep(0, p_tot*p_tot)
          for (k in 1:rep_num) {
            tmp <- co_glasso_D(corr_matrixes[[k]], p, t(as.matrix(c(alpha, coglasso_obj$lambda_w[j], lb_sel, c))), FALSE, FALSE, FALSE)
            convergence <- tmp$convergence[1]
            tmp <- as.vector(tmp$path[[1]])
            if (convergence == 1) {
              real_rep.num <- real_rep.num + 1
              merge_tmp <- merge_tmp + tmp
            }
          }
          rm(tmp)
          gc()
          
          merge_tmp <- merge_tmp / real_rep.num
          variability_tmp <- 4 * sum(merge_tmp * (1 - merge_tmp)) / (p_tot * (p_tot - 1))
          
          # Find way to deal with when real_rep.num is 0 (no convergence at all)
          if (variability_tmp >= stars_thresh) {
            index_lw <- max(j - 1, 1)
            tmp_lw <- coglasso_obj$lambda_w[[index_lw]]
            if (tmp_lw == lw_sel) converged <- TRUE
            
            lw_sel <- tmp_lw
            
            if (j == 1) {
              coglasso_obj$opt_variability[i] <- variability_tmp
              if (!light) {
                coglasso_obj$merge_lw[[i]] <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
              }
            }
            
            if (verbose) {
              mes <- "Conducting Subsampling Lambda_w....done.                 "
              cat(mes, "\r")
              cat("\n")
              if (converged) {
                mes <- "Reached convergence.                 "
                cat(mes, "\r")
                cat("\n")
              }
              flush.console()
            }
            
            break
          }
          
          else if (j == n_lambda_w) {
            index_lw <- 1
            tmp_lw <- coglasso_obj$lambda_w[[index_lw]]
            if (tmp_lw == lw_sel) converged <- TRUE
            
            lw_sel <- tmp_lw
            
            if (j == 1) {
              coglasso_obj$opt_variability[i] <- variability_tmp
              if (!light) {
                coglasso_obj$merge_lw[[i]] <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
              }
            }
            
            if (verbose) {
              mes <- "Conducting Subsampling Lambda_w....done.                 "
              cat(mes, "\r")
              cat("\n")
              if (converged) {
                mes <- "Reached convergence.                 "
                cat(mes, "\r")
                cat("\n")
              }
              flush.console()
            }
            
            break
          }
          
          coglasso_obj$opt_variability[i] <- variability_tmp
          if (!light) {
            coglasso_obj$merge_lw[[i]] <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
        }
      } 
      else {
        if (old_sampling) {
          corr_matrixes <- vector(mode = "list", length = rep_num)
          for (j in 1:rep_num) {
            ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
            corr_matrixes[[j]] <- cor(scale(coglasso_obj$data[ind.sample, ]))
          }
        }
        
        for (j in 1:n_lambda_b) {
          if (verbose) {
            mes <- paste(c("Conducting Subsampling Lambda_b....in progress:", floor(100 * j / n_lambda_b), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
          }
          real_rep.num <- 0
          merge_tmp <- rep(0, p_tot*p_tot)
          for (k in 1:rep_num) {
            tmp <- co_glasso_D(corr_matrixes[[k]], p, t(as.matrix(c(alpha, lw_sel, coglasso_obj$lambda_b[j], c))), FALSE, FALSE, FALSE)
            convergence <- tmp$convergence[1]
            tmp <- as.vector(tmp$path[[1]])
            if (convergence == 1) {
              real_rep.num <- real_rep.num + 1
              merge_tmp <- merge_tmp + tmp
            }
          }
          rm(tmp)
          gc()
          
          merge_tmp <- merge_tmp / real_rep.num
          variability_tmp <- 4 * sum(merge_tmp * (1 - merge_tmp)) / (p_tot * (p_tot - 1))
          
          if (variability_tmp >= stars_thresh) {
            index_lb <- max(j - 1, 1)
            tmp_lb <- coglasso_obj$lambda_b[[index_lb]]
            if (tmp_lb == lb_sel) converged <- TRUE
            
            lb_sel <- tmp_lb
            
            if (j == 1) {
              coglasso_obj$opt_variability[i] <- variability_tmp
              if (!light) {
                coglasso_obj$merge_lb[[i]] <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
              }
            }
            
            if (verbose) {
              mes <- "Conducting Subsampling Lambda_b....done.                 "
              cat(mes, "\r")
              cat("\n")
              if (converged) {
                mes <- "Reached convergence.                 "
                cat(mes, "\r")
                cat("\n")
              }
              flush.console()
            }
            
            break
          }
          
          else if (j == n_lambda_b) {
            index_lb <- 1
            tmp_lb <- coglasso_obj$lambda_b[[index_lb]]
            if (tmp_lb == lb_sel) converged <- TRUE
            
            lb_sel <- tmp_lb
            
            if (j == 1) {
              coglasso_obj$opt_variability[i] <- variability_tmp
              if (!light) {
                coglasso_obj$merge_lb[[i]] <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
              }
            }
            
            if (verbose) {
              mes <- "Conducting Subsampling Lambda_b....done.                 "
              cat(mes, "\r")
              cat("\n")
              if (converged) {
                mes <- "Reached convergence.                 "
                cat(mes, "\r")
                cat("\n")
              }
              flush.console()
            }
            
            break
          }
          
          coglasso_obj$opt_variability[i] <- variability_tmp
          if (!light) {
            coglasso_obj$merge_lb[[i]] <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
        }
      }
      
      select_lw <- !select_lw
      
      num_iter <- num_iter + 0.5
    }
    if (num_iter == max_iter) {
      if (verbose) {
        mes <- "Reached max Iterations.                 "
        cat(mes, "\r")
        cat("\n")
        flush.console()
      }
    }
    
    if (verbose) {
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values with \"xestars\"....done"
      cat(mes, "\r")
      flush.console()
    }
    
    coglasso_obj$opt_index_lw[i] <- index_lw
    coglasso_obj$opt_index_lb[i] <- index_lb
    coglasso_obj$opt_lambda_w[i] <- coglasso_obj$lambda_w[index_lw]
    coglasso_obj$opt_lambda_b[i] <- coglasso_obj$lambda_b[index_lb]
    opt_hpars_combination <- which(coglasso_obj$hpars[, 1] == alpha &
                                     coglasso_obj$hpars[, 2] == coglasso_obj$opt_lambda_w[i] &
                                     coglasso_obj$hpars[, 3] == coglasso_obj$opt_lambda_b[i] &
                                     coglasso_obj$hpars[, 4] == c)
    coglasso_obj$opt_adj[[i]] <- coglasso_obj$path[[opt_hpars_combination]]
    colnames(coglasso_obj$opt_adj[[i]]) <- colnames(coglasso_obj$data)
    row.names(coglasso_obj$opt_adj[[i]]) <- colnames(coglasso_obj$data)
  }
  
  if (all(is.na(coglasso_obj$opt_variability))) {
    warning("coglasso did not converge for any c parameter. Will select the highest c.")
  }
  
  coglasso_obj$sel_index_c <- min(which.min(coglasso_obj$opt_variability), n_c)
  coglasso_obj$sel_index_lw <- coglasso_obj$opt_index_lw[coglasso_obj$sel_index_c]
  coglasso_obj$sel_index_lb <- coglasso_obj$opt_index_lb[coglasso_obj$sel_index_c]
  coglasso_obj$sel_lambda_w <- coglasso_obj$lambda_w[coglasso_obj$sel_index_lw]
  coglasso_obj$sel_lambda_b <- coglasso_obj$lambda_b[coglasso_obj$sel_index_lb]
  coglasso_obj$sel_c <- coglasso_obj$c[coglasso_obj$sel_index_c]
  coglasso_obj$sel_adj <- coglasso_obj$opt_adj[[coglasso_obj$sel_index_c]]
  sel_hpars_combination <- which(coglasso_obj$hpars[, 1] ==  1 / (coglasso_obj$sel_c * (D - 1) + 1) &
                                   coglasso_obj$hpars[, 2] == coglasso_obj$sel_lambda_w &
                                   coglasso_obj$hpars[, 3] == coglasso_obj$sel_lambda_b &
                                   coglasso_obj$hpars[, 4] == coglasso_obj$sel_c)
  coglasso_obj$sel_density <- coglasso_obj$density[sel_hpars_combination]
  coglasso_obj$sel_icov <- coglasso_obj$icov[[sel_hpars_combination]]
  if (!is.null(coglasso_obj$cov)) {
    coglasso_obj$sel_cov <- coglasso_obj$cov[[sel_hpars_combination]]
  }
  coglasso_obj$method <- "xestars"
  coglasso_obj$call <- call
  
  class(coglasso_obj) <- "select_coglasso"
  
  return(coglasso_obj)
}

#' Stability selection of the best `coglasso` network
#' 
#' @description
#' `r lifecycle::badge("deprecated")`
#' 
#' `stars_coglasso()` was deprecated in favor of [xstars()] as the function name
#' does not reflect properly the method it implements. The new function will 
#' also allow new memory-lighter and faster possibilities. 
#' 
#' @keywords internal
#' @export
#' 
#' @examples
#' cg <- coglasso(multi_omics_sd_micro, p = c(4, 2), nlambda_w = 3, 
#'                nlambda_b = 3, nc = 3, verbose = FALSE)
#' \donttest{
#' # Deprecated, use xstars() instead. Takes around one minute
#' sel_cg <- stars_coglasso(cg, verbose = FALSE)
#' # ->
#' sel_cg <- xstars(cg, verbose = FALSE)
#' }
stars_coglasso <- function(coglasso_obj, stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, max_iter = 10, verbose = FALSE) {
  lifecycle::deprecate_warn("1.1.0", "stars_coglasso()", "xstars()", always = TRUE)
  xstars(coglasso_obj, stars_thresh, stars_subsample_ratio, rep_num, max_iter, verbose = FALSE)
}
