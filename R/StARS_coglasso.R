#' Stability selection of the best `coglasso` network
#'
#' `stars_coglasso()` selects the combination of hyperparameters given to
#' `coglasso()` yielding the most stable, yet sparse network. Stability is
#' computed upon network estimation from subsamples of the multi-omics data set,
#' allowing repetition. Subsamples are collected for a fixed amount of times
#' (`rep_num`), and with a fixed proportion of the total number of samples
#' (`stars_subsample_ratio`).
#'
#' *StARS* for *collaborative graphical regression* is an adaptation of the method
#' published by Liu, H. *et al.* (2010): Stability Approach to Regularization
#' Selection (StARS). *StARS* was developed for network estimation regulated by
#' a single penalty parameter, while collaborative graphical lasso needs to
#' explore three different hyperparameters. In particular, two of these are
#' penalty parameters with a direct influence on network sparsity, hence on
#' stability. For every \eqn{c} parameter, `stars_coglasso()` explores one of
#' the two penalty parameters (\eqn{\lambda_w} or \eqn{\lambda_b}), keeping the other one
#' fixed at its previous best estimate, using the normal, one-dimentional
#' *StARS* approach, until finding the best couple. It then selects the \eqn{c}
#' parameter for which the best (\eqn{\lambda_w}, \eqn{\lambda_b}) couple yielded the most
#' stable, yet sparse network.
#'
#' @param coglasso_obj The object returned by `coglasso()`.
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
#' @return `stars_coglasso()` returns a list containing the results of the
#'   selection procedure, built upon the list returned by `coglasso()`.
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
#'
#' @export
#'
#' @examples
#' cg <- coglasso(multi_omics_sd_micro, pX = 4, nlambda_w = 3, nlambda_b = 3, nc = 3, verbose = FALSE)
#' \dontrun{
#' # Takes around 20 seconds
#' sel_cg <- stars_coglasso(cg, verbose = FALSE)
#' }
stars_coglasso <- function(coglasso_obj, stars_thresh = 0.1, stars_subsample_ratio = NULL, rep_num = 20, max_iter = 10, verbose = TRUE) {
  n <- nrow(coglasso_obj$data)
  d <- ncol(coglasso_obj$data)
  pX <- coglasso_obj$pX
  n.lambda_w <- length(coglasso_obj$lambda_w)
  n.lambda_b <- length(coglasso_obj$lambda_b)
  n.c <- length(coglasso_obj$c)

  if (is.null(stars_subsample_ratio)) {
    if (n > 144) stars_subsample_ratio <- 10 * sqrt(n) / n
    if (n <= 144) stars_subsample_ratio <- 0.8
  }

  coglasso_obj$merge_lw <- vector(mode = "list", length = n.c)
  coglasso_obj$merge_lb <- vector(mode = "list", length = n.c)
  coglasso_obj$variability_lw <- vector(mode = "list", length = n.c)
  coglasso_obj$variability_lb <- vector(mode = "list", length = n.c)
  coglasso_obj$opt_adj <- vector(mode = "list", length = n.c)
  coglasso_obj$opt_variability <- rep(0, n.c)
  coglasso_obj$opt_index_lw <- rep(0, n.c)
  coglasso_obj$opt_index_lb <- rep(0, n.c)
  coglasso_obj$opt_lambda_w <- rep(0, n.c)
  coglasso_obj$opt_lambda_b <- rep(0, n.c)

  for (i in 1:n.c) {
    if (verbose) {
      mes <- paste(c("Selecting best Lambda_w/Lambda_b combination for all c values with \"stars\"....in progress:", floor(100 * i / n.c), "%"), collapse = "")
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }

    alpha <- 1 / (1 + coglasso_obj$c[i])

    lw_sel <- -1
    lb_sel <- coglasso_obj$lambda_b[n.lambda_b]
    num_iter <- 0
    converged <- FALSE
    select_lw <- TRUE

    while (!converged & num_iter < max_iter) {
      if (select_lw) {
        coglasso_obj$merge_lw[[i]] <- list()
        for (j in 1:n.lambda_w) coglasso_obj$merge_lw[[i]][[j]] <- Matrix::Matrix(0, d, d)

        real_rep.num <- rep(0, n.lambda_w)

        for (j in 1:rep_num)
        {
          if (verbose) {
            mes <- paste(c("Conducting Subsampling Lambda_w....in progress:", floor(100 * j / rep_num), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
          }
          ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)

          tmp <- coglasso(coglasso_obj$data[ind.sample, ], pX = pX, lambda_w = coglasso_obj$lambda_w, lambda_b = c(lb_sel), c = c(coglasso_obj$c[i]), verbose = FALSE)
          convergence <- tmp$convergence
          tmp <- tmp$path

          for (k in 1:n.lambda_w) {
            if (convergence[k] == 1) {
              real_rep.num[k] <- real_rep.num[k] + 1
              coglasso_obj$merge_lw[[i]][[k]] <- coglasso_obj$merge_lw[[i]][[k]] + tmp[[k]]
            }
          }
          rm(ind.sample, tmp)
          gc()
        }

        coglasso_obj$variability_lw[[i]] <- rep(0, n.lambda_w)
        for (j in 1:n.lambda_w) {
          coglasso_obj$merge_lw[[i]][[j]] <- coglasso_obj$merge_lw[[i]][[j]] / real_rep.num[j]
          coglasso_obj$variability_lw[[i]][j] <- 4 * sum(coglasso_obj$merge_lw[[i]][[j]] *
            (1 - coglasso_obj$merge_lw[[i]][[j]])) / (d * (d - 1))
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
        for (j in 1:n.lambda_b) coglasso_obj$merge_lb[[i]][[j]] <- Matrix::Matrix(0, d, d)

        real_rep.num <- rep(0, n.lambda_b)

        for (j in 1:rep_num)
        {
          if (verbose) {
            mes <- paste(c("Conducting Subsampling Lambda_b....in progress:", floor(100 * j / rep_num), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
          }
          ind.sample <- sample(c(1:n), floor(n * stars_subsample_ratio), replace = FALSE)

          tmp <- coglasso(coglasso_obj$data[ind.sample, ], pX = pX, lambda_w = c(lw_sel), lambda_b = coglasso_obj$lambda_b, c = c(coglasso_obj$c[i]), verbose = FALSE)
          convergence <- tmp$convergence
          tmp <- tmp$path

          for (k in 1:n.lambda_b) {
            if (convergence[k] == 1) {
              real_rep.num[k] <- real_rep.num[k] + 1
              coglasso_obj$merge_lb[[i]][[k]] <- coglasso_obj$merge_lb[[i]][[k]] + tmp[[k]]
            }
          }
          rm(ind.sample, tmp)
          gc()
        }

        coglasso_obj$variability_lb[[i]] <- rep(0, n.lambda_b)
        for (j in 1:n.lambda_b) {
          coglasso_obj$merge_lb[[i]][[j]] <- coglasso_obj$merge_lb[[i]][[j]] / real_rep.num[j]
          coglasso_obj$variability_lb[[i]][j] <- 4 * sum(coglasso_obj$merge_lb[[i]][[j]] *
            (1 - coglasso_obj$merge_lb[[i]][[j]])) / (d * (d - 1))
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
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values with \"stars\"....done"
      cat(mes, "\r")
      flush.console()
    }

    coglasso_obj$opt_index_lw[i] <- index_lw
    coglasso_obj$opt_index_lb[i] <- index_lb
    coglasso_obj$opt_lambda_w[i] <- coglasso_obj$lambda_w[index_lw]
    coglasso_obj$opt_lambda_b[i] <- coglasso_obj$lambda_b[index_lb]
    opt_hpars_combination <- which(coglasso_obj$hpars[, 1] == alpha &
      coglasso_obj$hpars[, 2] == coglasso_obj$opt_lambda_w[i] &
      coglasso_obj$hpars[, 3] == coglasso_obj$opt_lambda_b[i])
    coglasso_obj$opt_adj[[i]] <- coglasso_obj$path[[opt_hpars_combination]]
  }

  coglasso_obj$sel_index_c <- which.min(coglasso_obj$opt_variability)
  coglasso_obj$sel_index_lw <- coglasso_obj$opt_index_lw[coglasso_obj$sel_index_c]
  coglasso_obj$sel_index_lb <- coglasso_obj$opt_index_lb[coglasso_obj$sel_index_c]
  coglasso_obj$sel_lambda_w <- coglasso_obj$lambda_w[coglasso_obj$sel_index_lw]
  coglasso_obj$sel_lambda_b <- coglasso_obj$lambda_b[coglasso_obj$sel_index_lb]
  coglasso_obj$sel_c <- coglasso_obj$c[coglasso_obj$sel_index_c]
  coglasso_obj$sel_adj <- coglasso_obj$opt_adj[[coglasso_obj$sel_index_c]]
  sel_hpars_combination <- which(coglasso_obj$hpars[, 1] == 1 / (1 + coglasso_obj$sel_c) &
    coglasso_obj$hpars[, 2] == coglasso_obj$sel_lambda_w &
    coglasso_obj$hpars[, 3] == coglasso_obj$sel_lambda_b)
  coglasso_obj$sel_density <- coglasso_obj$density[sel_hpars_combination]
  coglasso_obj$sel_icov <- coglasso_obj$icov[[sel_hpars_combination]]
  if (!is.null(coglasso_obj$cov)) {
    coglasso_obj$sel_cov <- coglasso_obj$cov[[sel_hpars_combination]]
  }

  return(coglasso_obj)
}
