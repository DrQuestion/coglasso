#' Stability selection of the best `coglasso` network
#'
#' `stars_coglasso()` selects the combination of hyperparameters given to
#' `coglasso()` yielding the most stable, yet sparse network. Stability is
#' computed upon network estimation from subsamples of the multi-omics data set,
#' allowing repetition. Subsamples are collected for a fixed amount of times
#' (`rep.num`), and with a fixed proportion of the total number of samples
#' (`stars.subsample.ratio`).
#'
#' *StARS* for *collaborative graphical regression* is an adaptation of the method
#' published by Liu, H. *et al.* (2010): Stability Approach to Regularization
#' Selection (StARS). *StARS* was developed for network estimation regulated by
#' a single penalty parameter, while collaborative graphical lasso needs to
#' explore three different hyperparameters. In particular, two of these are
#' penalty parameters with a direct influence on network sparsity, hence on
#' stability. For every \eqn{c} parameter, `stars_coglasso()` explores one of
#' the two penalty parameters (\eqn{λ_w} or \eqn{λ_b}), keeping the other one
#' fixed at its previous best estimate, using the normal, one-dimentional
#' *StARS* approach, until finding the best couple. It then selects the \eqn{c}
#' parameter for which the best (\eqn{λ_w}, \eqn{λ_b}) couple yielded the most
#' stable, yet sparse network.
#'
#' @param coglasso.obj The object returned by `coglasso()`.
#' @param stars.thresh The threshold set for variability of the explored
#'   networks at each iteration of the algorithm. The \eqn{λ_w} or the \eqn{λ_b}
#'   associated to the most stable network before the threshold is overcome is
#'   selected.
#' @param stars.subsample.ratio The proportion of samples in the multi-omics
#'   data set to be randomly subsampled to estimate the variability of the
#'   network under the given hyperparameters setting. Defaults to 80% when the
#'   number of samples is smaller than 144, otherwise it defaults to
#'   \eqn{\frac{10}{n}\sqrt{n}}.
#' @param rep.num The amount of subsamples of the multi-omics data set used to
#'   estimate the variability of the network under the given hyperparameters
#'   setting. Defaults to 20.
#' @param max.iter The greatest number of times the algorithm is allowed to
#'   choose a new best \eqn{λ_w}. Defaults to 10.
#' @param verbose Print information regarding the progress of the selection
#'   procedure on the console.
#'
#' @return `stars_coglasso()` returns a list containing the results of the
#'   selection procedure, built upon the list returned by `coglasso()`.
#' * ... are the same elements returned by [coglasso()].
#' * `merge_lw` and `merge_lb` are lists with as many elements as the number of
#'   \eqn{c} parameters explored. Every element is in turn a list of as many
#'   matrices as the number of \eqn{λ_w} (or \eqn{λ_b}) values explored. Each
#'   matrix is the "merged" adjacency matrix, the average of all the adjacency
#'   matrices estimated  for those specific \eqn{c} and \eqn{λ_w} (or \eqn{λ_b})
#'   values across all the subsampling in the last path explored before
#'   convergence, the one when the final combination of \eqn{λ_w} and \eqn{λ_b}
#'   is selected for the given \eqn{c} value.
#' * `variability_lw` and `variability_lb` are lists with as many elements as 
#'  the number of \eqn{c} parameters explored. Every element is a numeric vector 
#'  of as many items as the number of \eqn{λ_w} (or \eqn{λ_b}) values explored. 
#'  Each item is the variability of the network estimated for those specific 
#'  \eqn{c} and \eqn{λ_w} (or \eqn{λ_b}) values in the last path explored before 
#'  convergence, the one when the final combination of \eqn{λ_w} and \eqn{λ_b} 
#'  is selected for the given \eqn{c} value.
#' * `opt.adj` is a list of the adjacency matrices finally selected for each 
#'  \eqn{c} parameter explored.
#' * `opt.variability` is a numerical vector containing the variabilities 
#'  associated to the adjacency matrices in `opt.adj`.
#' * `opt.index_lw` and `opt.index_lb` are integer vectors containing the 
#'  index of the selected \eqn{λ_w}s (or \eqn{λ_b}s) for each \eqn{c} parameters 
#'  explored.
#' * `opt.lambda_w` and `opt.lambda_b` are vectors containing the selected 
#'  \eqn{λ_w}s (or \eqn{λ_b}s) for each \eqn{c} parameters explored.
#' * `sel.index_c`, `sel.index_lw` and `sel.index_lb` are the indexes of the 
#'  final selected parameters \eqn{c}, \eqn{λ_w} and \eqn{λ_b} leading to the 
#'  most stable sparse network.
#' * `sel.c`, `sel.lambda_w` and `sel.lambda_b` are the final selected 
#'  parameters \eqn{c}, \eqn{λ_w} and \eqn{λ_b} leading to the most stable 
#'  sparse network.
#' * `sel.adj` is the adjacency matrix of the final selected network.
#' * `sel.density` is the density of the final selected network.
#' * `sel.icov` is the inverse covariance matrix of the final selected network.
#'
#' @export
#'
#' @examples
#' cg <- coglasso(multi_omics_SD_micro, pX=4, nlambda_w=6, nlambda_b=6, nc=3)
#' sel_cg <- stars_coglasso(cg)
#' 
stars_coglasso<-function(coglasso.obj, stars.thresh = 0.1, stars.subsample.ratio = NULL, rep.num = 20, max.iter = 10, verbose = TRUE)  
{
  n = nrow(coglasso.obj$data)
  d = ncol(coglasso.obj$data)
  pX <- coglasso.obj$pX
  n.lambda_w = length(coglasso.obj$lambda_w)
  n.lambda_b = length(coglasso.obj$lambda_b)
  n.c = length(coglasso.obj$c)
  
  if(is.null(stars.subsample.ratio))
  {
    if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
    if(n<=144) stars.subsample.ratio = 0.8
  }
  
  coglasso.obj$merge_lw = vector(mode = "list", length = n.c)
  coglasso.obj$merge_lb = vector(mode = "list", length = n.c)
  coglasso.obj$variability_lw = vector(mode = "list", length = n.c)
  coglasso.obj$variability_lb = vector(mode = "list", length = n.c)
  coglasso.obj$opt.adj = vector(mode = "list", length = n.c)
  coglasso.obj$opt.variability = rep(0, n.c)
  coglasso.obj$opt.index_lw = rep(0, n.c)
  coglasso.obj$opt.index_lb = rep(0, n.c)
  coglasso.obj$opt.lambda_w = rep(0, n.c)
  coglasso.obj$opt.lambda_b = rep(0, n.c)
  
  for (i in 1:n.c){
    
    if(verbose)
    {
      mes <- paste(c("Selecting best Lambda_w/Lambda_b combination for all c values....in progress:", floor(100*i/n.c), "%"), collapse="")
      cat(mes, "\r")
      cat("\n")
      flush.console()
    }
    
    alpha <- 1/(1+coglasso.obj$c[i])
    
    lw_sel <- -1
    lb_sel <- coglasso.obj$lambda_b[n.lambda_b]
    num_iter <- 0
    converged <- FALSE
    select_lw <- TRUE
    
    while (!converged & num_iter < max.iter) {
      if (select_lw){
        coglasso.obj$merge_lw[[i]] = list()
        for(j in 1:n.lambda_w) coglasso.obj$merge_lw[[i]][[j]] = Matrix::Matrix(0,d,d)
        
        real_rep.num <- rep(0, n.lambda_w)
        
        for(j in 1:rep.num)
        {
          if(verbose)
          {
            mes <- paste(c("Conducting Subsampling Lambda_w....in progress:", floor(100*j/rep.num), "%"), collapse="")
            cat(mes, "\r")
            flush.console()
          }
          ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)
          
          tmp = coglasso(coglasso.obj$data[ind.sample,], pX=pX, lambda_w = coglasso.obj$lambda_w, lambda_b = c(lb_sel), c=c(coglasso.obj$c[i]), verbose = FALSE)
          convergence <- tmp$convergence
          tmp <- tmp$path
          
          for(k in 1:n.lambda_w){
            if(convergence[k] == 1){
              real_rep.num[k] <- real_rep.num[k] + 1
              coglasso.obj$merge_lw[[i]][[k]] = coglasso.obj$merge_lw[[i]][[k]] + tmp[[k]]
            }
          }
          rm(ind.sample,tmp)
          gc()
        }
        
        coglasso.obj$variability_lw[[i]] = rep(0, n.lambda_w)
        for (j in 1:n.lambda_w) {
          coglasso.obj$merge_lw[[i]][[j]] = coglasso.obj$merge_lw[[i]][[j]]/real_rep.num[j]
          coglasso.obj$variability_lw[[i]][j] = 4 * sum(coglasso.obj$merge_lw[[i]][[j]] * 
                                                          (1 - coglasso.obj$merge_lw[[i]][[j]]))/(d * (d - 1))
        }
        
        index_lw = max(which.max(coglasso.obj$variability_lw[[i]] >= 
                                   stars.thresh)[1] - 1, 1)
        coglasso.obj$opt.variability[i] <- coglasso.obj$variability_lw[[i]][index_lw]
        tmp_lw <- coglasso.obj$lambda_w[[index_lw]]
        if(tmp_lw==lw_sel) converged <- TRUE
        
        lw_sel <- tmp_lw
        
        if (verbose) {
          mes = "Conducting Subsampling Lambda_w....done.                 "
          cat(mes, "\r")
          cat("\n")
          if (converged){
            mes = "Reached convergence.                 "
            cat(mes, "\r")
            cat("\n")
          }
          flush.console()
        }
      }
      
      else{
        coglasso.obj$merge_lb[[i]] = list()
        for(j in 1:n.lambda_b) coglasso.obj$merge_lb[[i]][[j]] = Matrix::Matrix(0,d,d)
        
        real_rep.num <- rep(0, n.lambda_b)
        
        for(j in 1:rep.num)
        {
          if(verbose)
          {
            mes <- paste(c("Conducting Subsampling Lambda_b....in progress:", floor(100*j/rep.num), "%"), collapse="")
            cat(mes, "\r")
            flush.console()
          }
          ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)
          
          tmp = coglasso(coglasso.obj$data[ind.sample,], pX=pX, lambda_w = c(lw_sel), lambda_b = coglasso.obj$lambda_b, c=c(coglasso.obj$c[i]), verbose = FALSE)
          convergence <- tmp$convergence
          tmp <- tmp$path
          
          for(k in 1:n.lambda_b){
            if(convergence[k] == 1){
              real_rep.num[k] <- real_rep.num[k] + 1
              coglasso.obj$merge_lb[[i]][[k]] = coglasso.obj$merge_lb[[i]][[k]] + tmp[[k]]
            }
            
          }
          rm(ind.sample,tmp)
          gc()
        }
        
        coglasso.obj$variability_lb[[i]] = rep(0, n.lambda_b)
        for (j in 1:n.lambda_b) {
          coglasso.obj$merge_lb[[i]][[j]] = coglasso.obj$merge_lb[[i]][[j]]/real_rep.num[j]
          coglasso.obj$variability_lb[[i]][j] = 4 * sum(coglasso.obj$merge_lb[[i]][[j]] * 
                                                          (1 - coglasso.obj$merge_lb[[i]][[j]]))/(d * (d - 1))
        }
        
        index_lb = max(which.max(coglasso.obj$variability_lb[[i]] >= 
                                   stars.thresh)[1] - 1, 1)
        coglasso.obj$opt.variability[i] <- coglasso.obj$variability_lb[[i]][index_lb]
        tmp_lb <- coglasso.obj$lambda_b[[index_lb]]
        if(tmp_lb==lb_sel) converged <- TRUE
        
        lb_sel <- tmp_lb
        
        if (verbose) {
          mes = "Conducting Subsampling Lambda_b....done.                 "
          cat(mes, "\r")
          cat("\n")
          if (converged){
            mes = "Reached convergence.                 "
            cat(mes, "\r")
            cat("\n")
          }
          flush.console()
        }
      }
      select_lw <- !select_lw
      
      num_iter <- num_iter + 0.5
    }
    if(num_iter==max.iter){
      if (verbose) {
        mes = "Reached max Iterations.                 "
        cat(mes, "\r")
        cat("\n")
        flush.console()
      }
    }
    
    
    if(verbose)
    {
      mes <- "Selecting best Lambda_w/Lambda_b combination for all c values....done"
      cat(mes, "\r")
      flush.console()
    }
    
    coglasso.obj$opt.index_lw[i] <- index_lw
    coglasso.obj$opt.index_lb[i] <- index_lb
    coglasso.obj$opt.lambda_w[i] <- coglasso.obj$lambda_w[index_lw]
    coglasso.obj$opt.lambda_b[i] <- coglasso.obj$lambda_b[index_lb]
    opt_hpars_combination <- which(coglasso.obj$hpars[,1]==alpha & 
                                     coglasso.obj$hpars[,2]==coglasso.obj$opt.lambda_w[i] & 
                                     coglasso.obj$hpars[,3]==coglasso.obj$opt.lambda_b[i])
    coglasso.obj$opt.adj[[i]] <- coglasso.obj$path[[opt_hpars_combination]]
  }
  
  coglasso.obj$sel.index_c <- which.min(coglasso.obj$opt.variability)
  coglasso.obj$sel.index_lw <- coglasso.obj$opt.index_lw[coglasso.obj$sel.index_c]
  coglasso.obj$sel.index_lb <- coglasso.obj$opt.index_lb[coglasso.obj$sel.index_c]
  coglasso.obj$sel.lambda_w <- coglasso.obj$lambda_w[coglasso.obj$sel.index_lw]
  coglasso.obj$sel.lambda_b <- coglasso.obj$lambda_b[coglasso.obj$sel.index_lb]
  coglasso.obj$sel.c <- coglasso.obj$c[coglasso.obj$sel.index_c]
  coglasso.obj$sel.adj <- coglasso.obj$opt.adj[[coglasso.obj$sel.index_c]]
  sel_hpars_combination <- which(coglasso.obj$hpars[,1]==1/(1+coglasso.obj$sel.c) & 
                                   coglasso.obj$hpars[,2]==coglasso.obj$sel.lambda_w & 
                                   coglasso.obj$hpars[,3]==coglasso.obj$sel.lambda_b)
  coglasso.obj$sel.density <- coglasso.obj$density[sel_hpars_combination]
  coglasso.obj$sel.icov = coglasso.obj$icov[[sel_hpars_combination]]
  if (!is.null(coglasso.obj$cov)) {
    coglasso.obj$sel.cov = coglasso.obj$cov[[sel_hpars_combination]]
  }
  
  return(coglasso.obj)
}
