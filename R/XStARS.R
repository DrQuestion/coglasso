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
#' explore three different hyperparameters. These all have, to a different 
#' degree, a direct influence on network sparsity, hence on stability. For every 
#' iteration, `xstars()` explores one of the three parameters (\eqn{\lambda_w}, 
#' \eqn{\lambda_b}, or \eqn{c}), keeping the other ones fixed at their previous 
#' selected estimate, using the normal, one-dimentional *StARS* approach, until 
#' finding the best combination of the three that yields the most stable, yet 
#' sparse network.
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
#' * `merge` is the "merged" adjacency matrix, the average of all the adjacency 
#'   matrices estimated across all the different subsamples for the selected 
#'   combination of \eqn{\lambda_w}, \eqn{\lambda_b}, and \eqn{c} values in the
#'   last path explored before convergence. Each entry is a measure of how 
#'   recurrent the corresponding edge is across the subsamples.
#' * `variability_lw`, `variability_lb` and `variability_c` are numeric vectors
#'   of as many items as the number of \eqn{\lambda_w}, \eqn{\lambda_b}, and 
#'   \eqn{c} values explored. Each item is the variability of the network 
#'   estimated for the corresponding hyperparameter value, keeping the other two 
#'   hyperparameters fixed to their selected value.
#' * `sel_index_c`, `sel_index_lw` and `sel_index_lb` are the indexes of the
#'   final selected parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} 
#'   leading to the most stable sparse network.
#' * `sel_c`, `sel_lambda_w` and `sel_lambda_b` are the final selected
#'   hyperparameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} leading to the 
#'   most stable sparse network.
#' * `sel_adj` is the adjacency matrix of the final selected network.
#' * `sel_variability` is the variability of the final selected network.
#' * `sel_density` is the density of the final selected network.
#' * `sel_icov` is the inverse covariance matrix of the final selected network.
#' * `sel_cov` optional, given only when `coglasso()` was called with 
#'   `cov_output = TRUE`. It is the covariance matrix associated with the final 
#'   selected network.
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
  
  coglasso_obj$variability_lw <- rep(0, n_lambda_w)
  coglasso_obj$variability_lb <- rep(0, n_lambda_b)
  coglasso_obj$variability_c <- rep(0, n_c)
  coglasso_obj$sel_adj <- Matrix::Matrix(0, p_tot, p_tot)
  coglasso_obj$sel_variability <- 0
  coglasso_obj$sel_index_lw <- 0
  coglasso_obj$sel_index_lb <- 0
  coglasso_obj$sel_index_c <- 0
  coglasso_obj$sel_lambda_w <- 0
  coglasso_obj$sel_lambda_b <- 0
  coglasso_obj$sel_c <- 0
  
  if (is.null(coglasso_obj$icov_guess)) {
    icov_guess <- matrix(0, p_tot, p_tot)
  }
  else {
    icov_guess <- coglasso_obj$icov_guess
  }
  
  if (verbose) {
    mes <- "Selecting best Lambda_w/Lambda_b/c combination with \"xstars\"....in progress"
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
  
  lw_sel <- -1
  lb_sel <- coglasso_obj$lambda_b[n_lambda_b]
  c_sel <- coglasso_obj$c[n_c]
  num_iter <- 0
  converged <- 0
  select_lw <- TRUE
  select_lb <- FALSE
  
  while ((converged < 2 | num_iter < 1) & num_iter < max_iter) {
    # XStARS converges when two selected hyperparameters in a row are equal to
    # the result of the previous selection procedure
    if (select_lw) {
      alpha <- 1 / (c_sel * (D - 1) + 1)
      
      coglasso_obj$merge <- vector(mode = "list", length = n_lambda_w)
      coglasso_obj$variability_lw <- rep(0, n_lambda_w)
      
      for (j in 1:n_lambda_w) coglasso_obj$merge[[j]] <- Matrix::Matrix(0, p_tot, p_tot)
      
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
        if (coglasso_obj$D == 2) {
          #hpars <- matrix(c(rep(alpha, n_lambda_w), coglasso_obj$lambda_w, rep(lb_sel, n_lambda_w)), nrow = n_lambda_w, ncol = 3)
          hpars <- matrix(c(rep(alpha, n_lambda_w), coglasso_obj$lambda_w, rep(lb_sel, n_lambda_w), rep(c_sel, n_lambda_w)), nrow = n_lambda_w, ncol = 4)
          #print(hpars)
          tmp <- co_glasso(corr_matrix, p[1], hpars, icov_guess, FALSE, FALSE, FALSE)
        }
        else {
          hpars <- matrix(c(rep(alpha, n_lambda_w), coglasso_obj$lambda_w, rep(lb_sel, n_lambda_w), rep(c_sel, n_lambda_w)), nrow = n_lambda_w, ncol = 4)
          tmp <- co_glasso_D(corr_matrix, p, hpars, FALSE, FALSE, FALSE)
        }
        
        convergence <- tmp$convergence
        tmp <- tmp$path
        
        for (k in 1:n_lambda_w) {
          if (convergence[k] == 1) {
            real_rep.num[k] <- real_rep.num[k] + 1
            coglasso_obj$merge[[k]] <- coglasso_obj$merge[[k]] + tmp[[k]]
          }
        }
        rm(ind.sample, tmp)
        gc()
      }
      
      for (j in 1:n_lambda_w) {
        coglasso_obj$merge[[j]] <- coglasso_obj$merge[[j]] / real_rep.num[j]
        coglasso_obj$variability_lw[j] <- 4 * sum(coglasso_obj$merge[[j]] *
                                                    (1 - coglasso_obj$merge[[j]])) / (p_tot * (p_tot - 1))
      }
      
      index_lw <- max(which.max(coglasso_obj$variability_lw >=
                                  stars_thresh)[1] - 1, 1)
      coglasso_obj$sel_variability <- coglasso_obj$variability_lw[index_lw]
      tmp_lw <- coglasso_obj$lambda_w[index_lw]
      if (tmp_lw == lw_sel) {
        converged <- converged + 1
      } else {
        converged <- 0
      }
      
      lw_sel <- tmp_lw
      
      if (verbose) {
        mes <- "Conducting Subsampling Lambda_w....done.                 "
        cat(mes, "\r")
        cat("\n")
        if (converged == 2) {
          mes <- "Reached convergence.                 "
          cat(mes, "\r")
          cat("\n")
        }
        flush.console()
      }
      select_lw <- FALSE
      select_lb <- TRUE
      
      if (converged == 2) {
        coglasso_obj$merge <- coglasso_obj$merge[[index_lw]]
      }
      
    } else if (select_lb) {
      alpha <- 1 / (c_sel * (D - 1) + 1)
      
      coglasso_obj$merge <- vector(mode = "list", length = n_lambda_b)
      coglasso_obj$variability_lb <- rep(0, n_lambda_b)
      
      for (j in 1:n_lambda_b) coglasso_obj$merge[[j]] <- Matrix::Matrix(0, p_tot, p_tot)
      
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
        if (coglasso_obj$D == 2) {
          #hpars <- matrix(c(rep(alpha, n_lambda_b), coglasso_obj$lambda_b, rep(lw_sel, n_lambda_b)), nrow = n_lambda_b, ncol = 3)
          hpars <- matrix(c(rep(alpha, n_lambda_b), rep(lw_sel, n_lambda_b), coglasso_obj$lambda_b, rep(c_sel, n_lambda_b)), nrow = n_lambda_b, ncol = 4)
          tmp <- co_glasso(corr_matrix, p[1], hpars, icov_guess, FALSE, FALSE, FALSE)
        }
        else {
          hpars <- matrix(c(rep(alpha, n_lambda_b), rep(lw_sel, n_lambda_b), coglasso_obj$lambda_b, rep(c_sel, n_lambda_b)), nrow = n_lambda_b, ncol = 4)
          tmp <- co_glasso_D(corr_matrix, p, hpars, FALSE, FALSE, FALSE)
        }
        
        convergence <- tmp$convergence
        tmp <- tmp$path
        
        for (k in 1:n_lambda_b) {
          if (convergence[k] == 1) {
            real_rep.num[k] <- real_rep.num[k] + 1
            coglasso_obj$merge[[k]] <- coglasso_obj$merge[[k]] + tmp[[k]]
          }
        }
        rm(ind.sample, tmp)
        gc()
      }
      
      coglasso_obj$variability_lb <- rep(0, n_lambda_b)
      for (j in 1:n_lambda_b) {
        coglasso_obj$merge[[j]] <- coglasso_obj$merge[[j]] / real_rep.num[j]
        coglasso_obj$variability_lb[j] <- 4 * sum(coglasso_obj$merge[[j]] *
                                                    (1 - coglasso_obj$merge[[j]])) / (p_tot * (p_tot - 1))
      }
      
      index_lb <- max(which.max(coglasso_obj$variability_lb >=
                                  stars_thresh)[1] - 1, 1)
      coglasso_obj$sel_variability <- coglasso_obj$variability_lb[index_lb]
      tmp_lb <- coglasso_obj$lambda_b[index_lb]
      if (tmp_lb == lb_sel) {
        converged <- converged + 1
      } else {
        converged <- 0
      }
      
      lb_sel <- tmp_lb
      
      if (verbose) {
        mes <- "Conducting Subsampling Lambda_b....done.                 "
        cat(mes, "\r")
        cat("\n")
        if (converged == 2) {
          mes <- "Reached convergence.                 "
          cat(mes, "\r")
          cat("\n")
        }
        flush.console()
      }
      
      if (converged == 2) {
        coglasso_obj$merge <- coglasso_obj$merge[[index_lb]]
      }
      
      select_lb <- FALSE
    } else {
      # Select c
      alpha <- 1 / (coglasso_obj$c * (D - 1) + 1)
      
      coglasso_obj$merge <- vector(mode = "list", length = n_c)
      coglasso_obj$variability_c <- rep(0, n_c)
      
      for (j in 1:n_c) coglasso_obj$merge[[j]] <- Matrix::Matrix(0, p_tot, p_tot)
      
      real_rep.num <- rep(0, n_c)
      
      for (j in 1:rep_num)
      {
        if (verbose) {
          mes <- paste(c("Conducting Subsampling c....in progress:", floor(100 * j / rep_num), "%"), collapse = "")
          cat(mes, "\r")
          flush.console()
        }
        ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
        
        corr_matrix <- cor(scale(coglasso_obj$data[ind.sample, ]))
        if (coglasso_obj$D == 2) {
          #hpars <- matrix(c(rep(alpha, n_c), coglasso_obj$c, rep(lw_sel, n_c)), nrow = n_c, ncol = 3)
          hpars <- matrix(c(alpha, rep(lw_sel, n_c), rep(lb_sel, n_c), coglasso_obj$c), nrow = n_c, ncol = 4)
          tmp <- co_glasso(corr_matrix, p[1], hpars, icov_guess, FALSE, FALSE, FALSE)
        }
        else {
          hpars <- matrix(c(alpha, rep(lw_sel, n_c), rep(lb_sel, n_c), coglasso_obj$c), nrow = n_c, ncol = 4)
          tmp <- co_glasso_D(corr_matrix, p, hpars, FALSE, FALSE, FALSE)
        }
        
        convergence <- tmp$convergence
        tmp <- tmp$path
        
        for (k in 1:n_c) {
          if (convergence[k] == 1) {
            real_rep.num[k] <- real_rep.num[k] + 1
            coglasso_obj$merge[[k]] <- coglasso_obj$merge[[k]] + tmp[[k]]
          }
        }
        rm(ind.sample, tmp)
        gc()
      }
      
      for (j in 1:n_c) {
        coglasso_obj$merge[[j]] <- coglasso_obj$merge[[j]] / real_rep.num[j]
        coglasso_obj$variability_c[j] <- 4 * sum(coglasso_obj$merge[[j]] *
                                                   (1 - coglasso_obj$merge[[j]])) / (p_tot * (p_tot - 1))
      }
      
      index_c <- max(which.max(coglasso_obj$variability_c >=
                                 stars_thresh)[1] - 1, 1)
      coglasso_obj$sel_variability <- coglasso_obj$variability_c[index_c]
      tmp_c <- coglasso_obj$c[index_c]
      if (tmp_c == c_sel) {
        converged <- converged + 1
      } else {
        converged <- 0
      }
      
      c_sel <- tmp_c
      
      if (verbose) {
        mes <- "Conducting Subsampling c....done.                 "
        cat(mes, "\r")
        cat("\n")
        if (converged == 2) {
          mes <- "Reached convergence.                 "
          cat(mes, "\r")
          cat("\n")
        }
        flush.console()
      }
      
      if (converged == 2) {
        coglasso_obj$merge <- coglasso_obj$merge[[index_c]]
      }
      
      select_lw <- TRUE
    }
    
    num_iter <- num_iter + 1/3
  }
  
  if (num_iter >= max_iter) {
    if (verbose) {
      mes <- "Reached max Iterations.                 "
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }
  }
  
  if (verbose) {
    mes <- "Selecting best Lambda_w/Lambda_b/c combination with \"xstars\"....done"
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
  
  coglasso_obj$sel_index_lw <- index_lw
  coglasso_obj$sel_index_lb <- index_lb
  coglasso_obj$sel_index_c <- index_c
  coglasso_obj$sel_lambda_w <- coglasso_obj$lambda_w[index_lw]
  coglasso_obj$sel_lambda_b <- coglasso_obj$lambda_b[index_lb]
  coglasso_obj$sel_c <- coglasso_obj$c[index_c]
  alpha <- 1 / (coglasso_obj$sel_c * (D - 1) + 1)
  sel_hpars_combination <- which(coglasso_obj$hpars[, 1] == alpha &
                                   coglasso_obj$hpars[, 2] == coglasso_obj$sel_lambda_w &
                                   coglasso_obj$hpars[, 3] == coglasso_obj$sel_lambda_b &
                                   coglasso_obj$hpars[, 4] == coglasso_obj$sel_c)
  coglasso_obj$sel_adj <- coglasso_obj$path[[sel_hpars_combination]]
  colnames(coglasso_obj$sel_adj) <- colnames(coglasso_obj$data)
  row.names(coglasso_obj$sel_adj) <- colnames(coglasso_obj$data)
  
  #if (is.na(coglasso_obj$sel_variability)) {
  #  warning("coglasso did not converge. Will select the highest c.")
  #}
  
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
#' explore three different hyperparameters. These all have, to different 
#' degree, a direct influence on network sparsity, hence on
#' stability. For every 
#' iteration, `xstars()` explores one of the three parameters (\eqn{\lambda_w}, 
#' \eqn{\lambda_b}, or \eqn{c}), keeping the other ones fixed at their previous 
#' selected estimate, using the normal, one-dimentional *StARS* approach, until 
#' finding the best combination of the three. What makes it more efficient
#' than `xstars()` is the different way that the stability check is implemented 
#' in the two algorithms. In `xstars()` (and even in the original *StARS*), the 
#' stability check is performed, for example, for every \eqn{\lambda_w} value 
#' (or \eqn{\lambda_b}, or \eqn{c}), until all values are explored, and then it 
#' when the algorithm selects the one yielding the most stable, yet sparse 
#' network, and only then switching to the selection of the following 
#' hyperparameter. In `xestars()`, the stability check becomes a *stopping criterion*. 
#' The moment that the stability threshold is passed, the value of the 
#' hyperparameter currently being selected is fixed, and the switch to the next 
#' one happens immediately, without exploring the whole landscape. This reduces 
#' sensibly the number of iterations before convergence to a final network. \cr
#' The original *XStARS* computes a new subsampling for every time the algorithm
#' switches from optimizing \eqn{\lambda_w}, \eqn{\lambda_b}, or \eqn{c}. This 
#' does not allow to compare the hyperparameters on an equal ground, and can 
#' slow the selection down with bigger data set or a larger hyperparameter 
#' space. To allow a similar subsampling to `xstars()`, the `old_sampling` 
#' parameter has been implemented. If set to TRUE, the subsampling is similar to 
#' the one `xstars()` would perform. Otherwise, the subsampling is performed at 
#' the beginning of the algorithm once and for all its iterations.
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
#' @param verbose Print information regarding the progress of the selection
#'   procedure on the console.
#'
#' @return `xestars()` returns an object of `S3` class `select_coglasso` 
#'   containing the results of the selection
#'   procedure, built upon the object of `S3` class `coglasso` returned by `coglasso()`.
#' * ... are the same elements returned by [coglasso()].
#' * `merge` is the "merged" adjacency matrix, the average of all the adjacency 
#'   matrices estimated across all the different subsamples for the selected 
#'   combination of \eqn{\lambda_w}, \eqn{\lambda_b}, and \eqn{c} values in the
#'   last path explored before convergence. Each entry is a measure of how 
#'   recurrent the corresponding edge is across the subsamples.
#' * `variability_lw`, `variability_lb` and `variability_c` are numeric vectors
#'   of as many items as the number of \eqn{\lambda_w}, \eqn{\lambda_b}, and 
#'   \eqn{c} values explored. Each item is the variability of the network 
#'   estimated for the corresponding hyperparameter value, keeping the other two 
#'   hyperparameters fixed to their selected value.
#' * `sel_index_c`, `sel_index_lw` and `sel_index_lb` are the indexes of the
#'   final selected parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b}
#'   leading to the most stable sparse network.
#' * `sel_c`, `sel_lambda_w` and `sel_lambda_b` are the final selected
#'   parameters \eqn{c}, \eqn{\lambda_w} and \eqn{\lambda_b} leading to the most
#'   stable sparse network.
#' * `sel_adj` is the adjacency matrix of the final selected network.
#' * `sel_variability` is the variability of the final selected network.
#' * `sel_density` is the density of the final selected network.
#' * `sel_icov` is the inverse covariance matrix of the final selected network.
#' * `sel_cov` optional, given only when `coglasso()` was called with 
#'   `cov_output = TRUE`. It is the covariance matrix associated with  the final 
#'   selected network.
#' * `call` is the matched call.
#' * `method` is the chosen model selection method. Here, it is "xestars".
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
xestars <- function(coglasso_obj, stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, max_iter = 10, old_sampling = FALSE, verbose = TRUE) {
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
  
  coglasso_obj$variability_lw <- rep(0, n_lambda_w)
  coglasso_obj$variability_lb <- rep(0, n_lambda_b)
  coglasso_obj$variability_c <- rep(0, n_c)
  coglasso_obj$sel_adj <- Matrix::Matrix(0, p_tot, p_tot)
  coglasso_obj$sel_variability <- 0
  coglasso_obj$sel_index_lw <- 0
  coglasso_obj$sel_index_lb <- 0
  coglasso_obj$sel_index_c <- 0
  coglasso_obj$sel_lambda_w <- 0
  coglasso_obj$sel_lambda_b <- 0
  coglasso_obj$sel_c <- 0
  
  # Create the same samples of indexes for the whole set of iterations,
  # cancelled when `old_sampling` is true
  if (!old_sampling) {
    corr_matrixes <- vector(mode = "list", length = rep_num)
    for (j in 1:rep_num) {
      ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
      corr_matrixes[[j]] <- cor(scale(coglasso_obj$data[ind.sample, ]))
    }
  }
  
  if (is.null(coglasso_obj$icov_guess)) {
    icov_guess <- matrix(0, p_tot, p_tot)
  }
  else {
    icov_guess <- coglasso_obj$icov_guess
  }
  
  if (verbose) {
    mes <- "Selecting best Lambda_w/Lambda_b/c combination with \"xestars\"....in progress"
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
  
  lw_sel <- coglasso_obj$lambda_w[n_lambda_w]
  lb_sel <- coglasso_obj$lambda_b[n_lambda_b]
  c_sel <- -1
  num_iter <- 0
  converged <- 0
  select_lw <- TRUE
  select_lb <- FALSE
  
  while ((converged < 2 | num_iter < 1) & num_iter < max_iter) {
    if (select_lw) {
      # Consider breaking in a one-dimensional sub-function returning a lambda estars
      alpha <- 1 / (c_sel * (D - 1) + 1)
      
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
          if (coglasso_obj$D == 2) {
            #tmp <- co_glasso(corr_matrixes[[k]], p[1], t(as.matrix(c(alpha, coglasso_obj$lambda_w[j], lb_sel))), FALSE, FALSE, FALSE)
            tmp <- co_glasso(corr_matrixes[[k]], p[1], t(as.matrix(c(alpha, coglasso_obj$lambda_w[j], lb_sel, c_sel))), icov_guess, FALSE, FALSE, FALSE)
          }
          else {
            tmp <- co_glasso_D(corr_matrixes[[k]], p, t(as.matrix(c(alpha, coglasso_obj$lambda_w[j], lb_sel, c_sel))), FALSE, FALSE, FALSE)
          }
          nexploded <- tmp$nexploded[1]
          tmp <- as.vector(tmp$path[[1]])
          if (nexploded == 0) {
            real_rep.num <- real_rep.num + 1
            merge_tmp <- merge_tmp + tmp
          }
        }
        rm(tmp)
        gc()
        
        merge_tmp <- merge_tmp / real_rep.num
        variability_tmp <- 4 * sum(merge_tmp * (1 - merge_tmp)) / (p_tot * (p_tot - 1))
        
        coglasso_obj$variability_lw[j] <- variability_tmp
        
        if (variability_tmp >= stars_thresh) {
          index_lw <- max(j - 1, 1)
          tmp_lw <- coglasso_obj$lambda_w[index_lw]
          if (tmp_lw == lw_sel) {
            converged <- converged + 1
          } else {
            converged <- 0
          }
          
          lw_sel <- tmp_lw
          
          if (j == 1) {
            coglasso_obj$sel_variability <- variability_tmp
            coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
          
          if (verbose) {
            mes <- "Conducting Subsampling Lambda_w....done.                 "
            cat(mes, "\r")
            cat("\n")
            if (converged == 2) {
              mes <- "Reached convergence.                 "
              cat(mes, "\r")
              cat("\n")
            }
            flush.console()
          }
          
          break
        } else if (j == n_lambda_w) {
          index_lw <- 1
          tmp_lw <- coglasso_obj$lambda_w[index_lw]
          if (tmp_lw == lw_sel) {
            converged <- converged + 1
          } else {
            converged <- 0
          }
          
          lw_sel <- tmp_lw
          
          if (j == 1) {
            coglasso_obj$sel_variability <- variability_tmp
            coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
          
          if (verbose) {
            mes <- "Conducting Subsampling Lambda_w....done.                 "
            cat(mes, "\r")
            cat("\n")
            if (converged == 2) {
              mes <- "Reached convergence.                 "
              cat(mes, "\r")
              cat("\n")
            }
            flush.console()
          }
          
          break
        }
        
        coglasso_obj$sel_variability <- variability_tmp
        coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
      }
      select_lw <- FALSE
      select_lb <- TRUE
      
    } else if (select_lb) {
      alpha <- 1 / (c_sel * (D - 1) + 1)
      
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
          if (coglasso_obj$D == 2) {
            #tmp <- co_glasso(corr_matrixes[[k]], p[1], t(as.matrix(c(alpha, lw_sel, coglasso_obj$lambda_b[j]))), FALSE, FALSE, FALSE)
            tmp <- co_glasso(corr_matrixes[[k]], p[1], t(as.matrix(c(alpha, lw_sel, coglasso_obj$lambda_b[j], c_sel))), icov_guess, FALSE, FALSE, FALSE)
          }
          else {
            tmp <- co_glasso_D(corr_matrixes[[k]], p, t(as.matrix(c(alpha, lw_sel, coglasso_obj$lambda_b[j], c_sel))), FALSE, FALSE, FALSE)
          }
          nexploded <- tmp$nexploded[1]
          tmp <- as.vector(tmp$path[[1]])
          if (nexploded == 0) {
            real_rep.num <- real_rep.num + 1
            merge_tmp <- merge_tmp + tmp
          }
        }
        rm(tmp)
        gc()
        
        merge_tmp <- merge_tmp / real_rep.num
        variability_tmp <- 4 * sum(merge_tmp * (1 - merge_tmp)) / (p_tot * (p_tot - 1))
        
        coglasso_obj$variability_lb[j] <- variability_tmp
        
        if (variability_tmp >= stars_thresh) {
          index_lb <- max(j - 1, 1)
          tmp_lb <- coglasso_obj$lambda_b[index_lb]
          if (tmp_lb == lb_sel) {
            converged <- converged + 1
          } else {
            converged <- 0
          }
          
          lb_sel <- tmp_lb
          
          if (j == 1) {
            coglasso_obj$sel_variability <- variability_tmp
            coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
          
          if (verbose) {
            mes <- "Conducting Subsampling Lambda_b....done.                 "
            cat(mes, "\r")
            cat("\n")
            if (converged == 2) {
              mes <- "Reached convergence.                 "
              cat(mes, "\r")
              cat("\n")
            }
            flush.console()
          }
          
          break
        } else if (j == n_lambda_b) {
          index_lb <- 1
          tmp_lb <- coglasso_obj$lambda_b[index_lb]
          if (tmp_lb == lb_sel) {
            converged <- converged + 1
          } else {
            converged <- 0
          }
          
          lb_sel <- tmp_lb
          
          if (j == 1) {
            coglasso_obj$sel_variability <- variability_tmp
            coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
          
          if (verbose) {
            mes <- "Conducting Subsampling Lambda_b....done.                 "
            cat(mes, "\r")
            cat("\n")
            if (converged == 2) {
              mes <- "Reached convergence.                 "
              cat(mes, "\r")
              cat("\n")
            }
            flush.console()
          }
          
          break
        }
        
        coglasso_obj$sel_variability <- variability_tmp
        coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
      }
      select_lb <- FALSE
    } else  {
      if (old_sampling) {
        corr_matrixes <- vector(mode = "list", length = rep_num)
        for (j in 1:rep_num) {
          ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)
          corr_matrixes[[j]] <- cor(scale(coglasso_obj$data[ind.sample, ]))
        }
      }
      
      for (j in 1:n_c) {
        alpha <- 1 / (coglasso_obj$c[j] * (D - 1) + 1)
        
        if (verbose) {
          mes <- paste(c("Conducting Subsampling c....in progress:", floor(100 * j / n_c), "%"), collapse = "")
          cat(mes, "\r")
          flush.console()
        }
        real_rep.num <- 0
        merge_tmp <- rep(0, p_tot*p_tot)
        for (k in 1:rep_num) {
          if (coglasso_obj$D == 2) {
            #tmp <- co_glasso(corr_matrixes[[k]], p[1], t(as.matrix(c(alpha, lw_sel, coglasso_obj$lambda_b[j]))), FALSE, FALSE, FALSE)
            tmp <- co_glasso(corr_matrixes[[k]], p[1], t(as.matrix(c(alpha, lw_sel, lb_sel, coglasso_obj$c[j]))), icov_guess, FALSE, FALSE, FALSE)
          }
          else {
            tmp <- co_glasso_D(corr_matrixes[[k]], p, t(as.matrix(c(alpha, lw_sel, lb_sel, coglasso_obj$c[j]))), FALSE, FALSE, FALSE)
          }
          nexploded <- tmp$nexploded[1]
          tmp <- as.vector(tmp$path[[1]])
          if (nexploded == 0) {
            real_rep.num <- real_rep.num + 1
            merge_tmp <- merge_tmp + tmp
          }
        }
        rm(tmp)
        gc()
        
        merge_tmp <- merge_tmp / real_rep.num
        variability_tmp <- 4 * sum(merge_tmp * (1 - merge_tmp)) / (p_tot * (p_tot - 1))
        
        coglasso_obj$variability_c[j] <- variability_tmp
        
        if (variability_tmp >= stars_thresh) {
          index_c <- max(j - 1, 1)
          tmp_c <- coglasso_obj$c[index_c]
          if (tmp_c == c_sel) {
            converged <- converged + 1
          } else {
            converged <- 0
          }
          
          c_sel <- tmp_c
          
          if (j == 1) {
            coglasso_obj$sel_variability <- variability_tmp
            coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
          
          if (verbose) {
            mes <- "Conducting Subsampling c....done.                 "
            cat(mes, "\r")
            cat("\n")
            if (converged == 2) {
              mes <- "Reached convergence.                 "
              cat(mes, "\r")
              cat("\n")
            }
            flush.console()
          }
          
          break
        } else if (j == n_c) {
          index_c <- 1
          tmp_c <- coglasso_obj$c[index_c]
          if (tmp_c == c_sel) {
            converged <- converged + 1
          } else {
            converged <- 0
          }
          
          c_sel <- tmp_c
          
          if (j == 1) {
            coglasso_obj$sel_variability <- variability_tmp
            coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
          }
          
          if (verbose) {
            mes <- "Conducting Subsampling c....done.                 "
            cat(mes, "\r")
            cat("\n")
            if (converged == 2) {
              mes <- "Reached convergence.                 "
              cat(mes, "\r")
              cat("\n")
              cat("No c value was enough to pass the threshold of instability.\nTry with a higher value.\n")
              
            }
            flush.console()
          }
          
          break
        }
        
        coglasso_obj$sel_variability <- variability_tmp
        coglasso_obj$merge <- merge_tmp #make sure to make it a matrix again if implementing the vector to compute variability above
      }
      select_lw <- TRUE
    }
    
    num_iter <- num_iter + 1/3
  }
  if (num_iter >= max_iter) {
    if (verbose) {
      mes <- "Reached max Iterations.                 "
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }
  }
  
  if (verbose) {
    mes <- "Selecting best Lambda_w/Lambda_b/c combination with \"xestars\"....done"
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
  
  coglasso_obj$sel_index_lw <- index_lw
  coglasso_obj$sel_index_lb <- index_lb
  coglasso_obj$sel_index_c <- index_c
  coglasso_obj$sel_lambda_w <- coglasso_obj$lambda_w[index_lw]
  coglasso_obj$sel_lambda_b <- coglasso_obj$lambda_b[index_lb]
  coglasso_obj$sel_c <- coglasso_obj$c[index_c]
  alpha <- 1 / (coglasso_obj$sel_c * (D - 1) + 1)
  sel_hpars_combination <- which(coglasso_obj$hpars[, 1] == alpha &
                                   coglasso_obj$hpars[, 2] == coglasso_obj$sel_lambda_w &
                                   coglasso_obj$hpars[, 3] == coglasso_obj$sel_lambda_b &
                                   coglasso_obj$hpars[, 4] == coglasso_obj$sel_c)
  coglasso_obj$sel_adj <- coglasso_obj$path[[sel_hpars_combination]]
  colnames(coglasso_obj$sel_adj) <- colnames(coglasso_obj$data)
  row.names(coglasso_obj$sel_adj) <- colnames(coglasso_obj$data)
  
  #if (all(is.na(coglasso_obj$sel_variability))) {
  #  warning("coglasso did not converge for any c parameter. Will select the highest c.")
  #}
  
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
#' does not reflect properly the method it implements. The new function 
#' also allows new memory-lighter and faster possibilities. 
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
