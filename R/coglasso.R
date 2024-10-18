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
#' @param p A vector with with the number of variables for each omic layer of the 
#'   data set (e.g. the number of transcripts, metabolites etc.), in the same 
#'   order the layers have in the data set. If given a single number, 
#'   `coglasso()` assumes that the total of data sets is two, and that the 
#'   number given is the dimension of the first one.
#' @param pX `r lifecycle::badge("deprecated")` `pX` is no longer supported. 
#'   Please use `p`.
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
#' @param icov_guess Use a predetermined inverse covariance matrix as an initial
#'   guess for the network estimation.
#' @param cov_output Add the estimated variance-covariance matrix to the output.
#' @param lock_lambdas Set \eqn{\lambda_w = \lambda_b}. Force a single lambda 
#'   parameter for both "within" and "between" interactions.
#' @param verbose Print information regarding current `coglasso` run on the
#'   console.
#'
#' @return `coglasso()` returns an object of `S3` class `coglasso`, that has the 
#' following elements:
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
#' * `p` is the vector with the number of variables for each omic layer of the 
#'   data set.
#' * `D` is the number of omics layers in the data set.
#' * `icov_guess` optional, returned when `icov_guess` is given. It is the 
#'   predetermined inverse covariance matrix  given by the user as an initial 
#'   guess for the network estimation.
#' * `cov` optional, returned when `cov_output` is TRUE, is a list containing
#'   the variance-covariance matrices of all the estimated networks.
#' * `call` is the matched call.
#'
#' @export
#'
#' @examples
#' # Typical usage: set the number of hyperparameters to explore
#' cg <- coglasso(multi_omics_sd_micro, p = c(4, 2), nlambda_w = 3, 
#'                nlambda_b = 3, nc = 3, verbose = FALSE)
#' \donttest{
#' # Model selection using eXtended Efficient StARS, takes less than five seconds
#' sel_cg_xestars <- select_coglasso(cg, method = "xestars", verbose = FALSE)
#' }
#' 
coglasso <- function(data, p = NULL, pX = lifecycle::deprecated(), lambda_w = NULL, lambda_b = NULL, c = NULL, nlambda_w = NULL, nlambda_b = NULL, nc = NULL, lambda_w_max = NULL, lambda_b_max = NULL, c_max = NULL, lambda_w_min_ratio = NULL, lambda_b_min_ratio = NULL, c_min_ratio = NULL, icov_guess = NULL, cov_output = FALSE, lock_lambdas = FALSE, verbose = TRUE) {
  call <- match.call()
  
  if (lifecycle::is_present(pX)) {
    lifecycle::deprecate_warn("1.1.0", "coglasso(pX)", "coglasso(p)")
    p <- pX
  }
  
  original_data <- data
  data <- as.matrix(data)
  gcinfo(FALSE)
  n <- nrow(data)
  p_tot <- ncol(data)
  
  if (is.null(p) && is.null(pX)) {
    stop("Please specify the dimension of the data sets.")
  }
  if (length(p)==1) {
    p <- c(p, p_tot - p)
  }
  
  if (is.null(icov_guess)) {
    icov_guess_null <- TRUE
    icov_guess <- matrix(0, p_tot, p_tot)
  }
  
  D <- length(p)
  
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
    S = S, p = p, D = D, lambda_w = lambda_w, lambda_b = lambda_b, c = c, nlambda_w = nlambda_w, nlambda_b = nlambda_b, nc = nc, lambda_w_max = lambda_w_max, lambda_b_max = lambda_b_max, c_max = c_max,
    lambda_w_min_ratio = lambda_w_min_ratio, lambda_b_min_ratio = lambda_b_min_ratio, c_min_ratio = c_min_ratio, lock_lambdas = lock_lambdas 
  )
  
  if (D == 2) {
    cg <- co_glasso(S, p[1], matrix(hpars[[1]], ncol = 4), icov_guess, FALSE, verbose, cov_output)
  }
  else {
    cg <- co_glasso_D(S, p, matrix(hpars[[1]], ncol = 4), FALSE, verbose, cov_output)
  }
  
  if (!icov_guess_null) {
    cg$icov_guess <- icov_guess
  }
  cg$data <- original_data
  cg$hpars <- hpars[[1]]
  cg$lambda_w <- hpars[[2]]
  cg$lambda_b <- hpars[[3]]
  cg$c <- hpars[[4]]
  # Modify for |D|:
  cg$p <- p
  cg$D <- D
  cg$call <- call
  
  class(cg) <- "coglasso"
  
  cg
}

#' Print function for the S3 class `coglasso`
#'
#' Print information on the estimated networks and the explored hyperparameters 
#' and see next suggested step
#'
#' @param x The object of `S3` class `coglasso`.
#' @param ... System required, not used.
#' 
#' @noRd
#' @export
print.coglasso <- function(x, ...){
  x_name <- rlang::call_args(match.call())
  cat("Networks estimated with collaborative graphical lasso\n\n")
  cat("The call was:\n", 
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("The total number of hyperparameter combinations explored was:\n", 
      dim(x$hpars)[1], "\n", sep = "")
  cat("The values explored for lambda within were:\n", 
      paste(round(x$lambda_w, 4), collapse = ", "), "\n", sep = "")
  cat("The values explored for lambda between were:\n", 
      paste(round(x$lambda_b, 4), collapse = ", "), "\n", sep = "")
  cat("The values explored for c were:\n", paste(round(x$c, 4), collapse = ", "), "\n\n",
      sep = "")
  cat("Networks are made of ", x$D, " omics layers, for a total of ", ncol(x$data), " nodes\n", sep = "")
  mes <- paste(x$p[-length(x$p)], collapse = ", ")
  cat("For each layer they have: ", mes, " and ", x$p[length(x$p)], " nodes, respectively\n\n", sep = "")
  cat("Select the best network with:\nsel_cg <- select_coglasso(", x_name[[1]], ")",
      sep = "")
}

#' Generate combinations of hyperparameters for `coglasso()`
#'
#' `gen_hpars()` generates an ordered table of all the combinations of the hyperparameters given as input to `coglasso()`
#'
#' @inherit coglasso
#' @param S The empirical Pearson's correlation matrix of the multi-omics data set.
#' @param D The number of omics layers in the data set.
#'
#' @return
#' `gen_hpars()` returns a list of four elements:
#' * A table of all the combinations of hyperparameters given as input to `coglasso()`, with \eqn{\alpha(\lambda_w+\lambda_b)} being the key to sort rows.
#' * A numerical vector with all the generated \eqn{\lambda_w}.
#' * A numerical vector with all the generated \eqn{\lambda_b}.
#' * A numerical vector with all the generated \eqn{c}.
#' 
#' @noRd
gen_hpars <- function(S = NULL, p = NULL, D = NULL, lambda_w = NULL, lambda_b = NULL, c = NULL,
                      nlambda_w = NULL, nlambda_b = NULL, nc = NULL, 
                      lambda_w_max = NULL, lambda_b_max = NULL, c_max = NULL, 
                      lambda_w_min_ratio = NULL, lambda_b_min_ratio = NULL,
                      c_min_ratio = NULL, lock_lambdas = FALSE) {
  # Modify alpha generation based on |D|
  # Add a fourth column containing c --> delete getc as an option
  if (!is.null(lambda_w)) nlambda_w <- length(lambda_w)
  if (!is.null(lambda_b)) nlambda_b <- length(lambda_b)
  if (!is.null(c)) nc <- length(c)
  
  p_tot <- sum(p)
  
  if (lock_lambdas) {
    condition1 <- !is.null(lambda_w) & !is.null(lambda_b)
    condition2 <- !is.null(nlambda_w) & !is.null(nlambda_b) 
    condition3 <- !is.null(lambda_w) & !is.null(nlambda_b)
    condition4 <- !is.null(nlambda_w) & !is.null(lambda_b)
    
    if (condition1) {
      if(lambda_w != lambda_b) {
        stop("Cannot lock lambda if different values are given for lambda_w and 
             lambda_b. Please either provide values for only one of the two 
             arguments or provide the same values")
      } 
    }
    else if (condition2 | condition3 | condition4) {
      if (nlambda_w != nlambda_b) {
        stop("Cannot lock lambda if the amount of values to explore for the two 
             lambdas differs. Please either provide values for only one of the 
             two arguments or provide coherent values")
      }
    }
    
    if (!is.null(lambda_w)) {
      lambda_b <- lambda_w
    }
    else if (!is.null(lambda_b)) {
      lambda_w <- lambda_b
    }
    else if (!is.null(nlambda_w)) {
      nlambda_b <- nlambda_w
    }
    else if (!is.null(nlambda_b)) {
      nlambda_w <- nlambda_b
    }
  }
  
  p_i <- 1
  S_diag_blocks <- matrix(0, p_tot, p_tot)
  for (i in 1:length(p)) {
    coords <- seq(p_i, p[i] + p_i - 1)
    S_diag_blocks[coords, coords] <- S[coords, coords]
    p_i <- p_i + p[i]
  }
  
  if (is.null(lambda_w)) {
    if (is.null(nlambda_w)) {
      nlambda_w <- 8
    }
    if (is.null(lambda_w_min_ratio)) {
      lambda_w_min_ratio <- 0.1
    }
    if (is.null(lambda_w_max)) {
      lambda_w_max <- max(max(S_diag_blocks - diag(p_tot)), 
                          -min(S_diag_blocks - diag(p_tot)))
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
      lambda_b_max <- max(max(S - S_diag_blocks), -min(S - S_diag_blocks))
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
  
  alpha <- 1 / (c * (D - 1) +1)
  ac <- cbind(alpha, c)
  if (lock_lambdas) {
    hpars_table <- expand.grid(alpha, lambda_w)
    hpars_table <- cbind(hpars_table, hpars_table[, 2])
  }
  else {
    hpars_table <- expand.grid(alpha, lambda_w, lambda_b)
  }
  hpars_table$id <- 1:nrow(hpars_table)
  hpars_table <- merge(hpars_table, ac, by = 1)
  hpars_table <- hpars_table[order(hpars_table$id), ]
  key <- hpars_table[, 1] * (hpars_table[, 2] + hpars_table[, 3])
  hpars_table <- cbind(hpars_table, key)
  hpars_table <- hpars_table[order(hpars_table$key, decreasing = T), ]
  hpars_table <- hpars_table[, c(1, 2, 3, 5)]
  
  hpars[[1]] <- as.matrix(hpars_table)
  colnames(hpars[[1]]) <- c("alpha", "lambda_w", "lambda_b", "c")
  rownames(hpars[[1]]) <- seq(1:nrow(hpars[[1]]))
  hpars[[2]] <- lambda_w
  hpars[[3]] <- lambda_b
  hpars[[4]] <- c
  
  hpars
}
