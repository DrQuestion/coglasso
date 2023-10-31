#' Estimate networks from a multi-omics data set
#'
#' `coglasso()` estimates multiple multi-omics networks with the algorithm
#' *Collaborative Graphical Lasso*, one for each combination of input values for
#' the hyperparameters \eqn{λ_w}, \eqn{λ_b} and \eqn{c}.
#'
#' @param data The input multi-omics data set. Rows should be samples, columns
#'   should be variables. Variables should be grouped by their assay (i.e.
#'   transcripts first, then metabolites). `data` is a required parameter.
#' @param pX The number of variables of the first data set (i.e. the number of
#'   transcripts). `pX` is a required parameter.
#' @param lambda_w A vector of values for the parameter \eqn{λ_w}, the
#'   penalization parameter for the "within" interactions. Overrides
#'   `nlambda_w`.
#' @param lambda_b A vector of values for the parameter \eqn{λ_b}, the
#'   penalization parameter for the "between" interactions. Overrides
#'   `nlambda_b`.
#' @param c A vector of values for the parameter \eqn{c}, the weight given to
#'   collaboration. Overrides `nc`.
#' @param nlambda_w The number of requested \eqn{λ_w} parameters to explore. A
#'   sequence of size `nlambda_w` of \eqn{λ_w} parameters will be generated.
#'   Defaults to 8. Ignored when `lambda_w` is set by the user.
#' @param nlambda_b The number of requested \eqn{λ_b} parameters to explore. A
#'   sequence of size `nlambda_b` of \eqn{λ_b} parameters will be generated.
#'   Defaults to 8. Ignored when `lambda_b` is set by the user.
#' @param nc The number of requested \eqn{c} parameters to explore. A sequence
#'   of size `nc` of \eqn{c} parameters will be generated. Defaults to 8.
#'   Ignored when `c` is set by the user.
#' @param c.max The greatest generated \eqn{c}. Defaults to 10. Ignored when `c`
#'   is set by the user.
#' @param lambda_w.min.ratio The ratio of the smallest generated \eqn{λ_w} over
#'   the greatest generated \eqn{λ_w}. Defaults to 0.1. Ignored when `lambda_w`
#'   is set by the user.
#' @param lambda_b.min.ratio The ratio of the smallest generated \eqn{λ_b} over
#'   the greatest generated \eqn{λ_b}. Defaults to 0.1. Ignored when `lambda_b`
#'   is set by the user.
#' @param c.min.ratio The ratio of the smallest generated \eqn{c} over the
#'   greatest generated \eqn{c}. Defaults to 0.1. Ignored when `c` is set by the
#'   user.
#' @param cov.output Add the estimated variance-covariance matrix to the output.
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
#'   given as input to `coglasso()`, with \eqn{α(λ_w+λ_b)} being the key to sort 
#'   rows.
#' * `lambda_w` is a numerical vector with all the \eqn{λ_w} values `coglasso()` 
#'   used.
#' * `lambda_b` is a numerical vector with all the \eqn{λ_b} values `coglasso()` 
#'   used.
#' * `c` is a numerical vector with all the \eqn{c} values `coglasso()` 
#'   used.
#' * `pX` is the number of variables of the first data set.
#' * `cov` optional, returned when `cov.output` is TRUE, is a list containing 
#'   the variance-covariance matrices of all the estimated networks.
#'
#' @export
#' 
#' @examples
#' # Typical usage: set the number of hyperparameters to explore
#' cg <- coglasso(multi_omics_SD_micro, pX=4, nlambda_w=6, nlambda_b=6, nc=3)
#' 
coglasso <- function(data, pX, lambda_w = NULL, lambda_b = NULL, c = NULL, nlambda_w = NULL, nlambda_b = NULL, nc = NULL, c.max=NULL, lambda_w.min.ratio = NULL, lambda_b.min.ratio = NULL, c.min.ratio = NULL, cov.output = FALSE, verbose = TRUE){
  original_data = data
  data = as.matrix(data)
  gcinfo(FALSE)
  n = nrow(data)
  p = ncol(data)
  cov.input = isSymmetric(data)
  if(cov.input)
  {
    if(verbose) cat("The input is identified as the covariance matrix.\n")
    S = data
  }
  else
  {
    data = scale(data)
    S = cor(data)
  }
  rm(data)
  gc()
  
  hpars<-gen_hpars(S=S, p=p, lambda_w=lambda_w, lambda_b=lambda_b, c=c, nlambda_w=nlambda_w, nlambda_b=nlambda_b, nc=nc, c.max=c.max, 
                   lambda_w.min.ratio=lambda_w.min.ratio, lambda_b.min.ratio=lambda_b.min.ratio, c.min.ratio=c.min.ratio)
  
  cg<-co_glasso(S, pX, hpars[[1]], FALSE, verbose, cov.output)
  
  cg$data <- original_data
  cg$hpars <- hpars[[1]]
  cg$lambda_w <- hpars[[2]]
  cg$lambda_b <- hpars[[3]]
  cg$c <- hpars[[4]]
  cg$pX <- pX
  
  cg
}

#' Generate combinations of hyperparameters for `coglasso()`
#' 
#' `gen_hpars()` generates an ordered table of all the combinations of the hyperparameters given as input to `coglasso()`
#'
#' @inherit coglasso
#' @param S The empirical Pearson's correlation matrix of the multi-omics data set.
#' @param p The total number of variables in the multi-omics data set.
#' @param getc Parameter for testing purposes. If TRUE, the returned table contains a column with *c* values rather than *α* values.
#'
#' @return
#' `gen_hpars()` returns a list of four elements: 
#' * A table of all the combinations of hyperparameters given as input to `coglasso()`, with \eqn{α(λ_w+λ_b)} being the key to sort rows.
#' * A numerical vector with all the generated \eqn{λ_w}.
#' * A numerical vector with all the generated \eqn{λ_b}.
#' * A numerical vector with all the generated \eqn{c}.
#' 
gen_hpars <- function(S=NULL, p=NULL, lambda_w = NULL, lambda_b = NULL, c = NULL, 
                      nlambda_w = NULL, nlambda_b = NULL, 
                      nc = NULL, c.max=NULL, lambda_w.min.ratio = NULL, lambda_b.min.ratio = NULL, 
                      c.min.ratio = NULL, getc = NULL){
  if(!is.null(lambda_w)) nlambda_w = length(lambda_w)
  if(!is.null(lambda_b)) nlambda_b = length(lambda_b)
  if(!is.null(c)) nalpha = length(c)
  
  if(is.null(lambda_w))
  {
    if(is.null(nlambda_w))
      nlambda_w = 8
    if(is.null(lambda_w.min.ratio))
      lambda_w.min.ratio = 0.1
    lambda_w.max = max(max(S-diag(p)),-min(S-diag(p)))
    lambda_w.min = lambda_w.min.ratio*lambda_w.max
    lambda_w = exp(seq(log(lambda_w.max), log(lambda_w.min), length = nlambda_w))
  }
  if(is.null(lambda_b))
  {
    if(is.null(nlambda_b))
      nlambda_b = 8
    if(is.null(lambda_b.min.ratio))
      lambda_b.min.ratio = 0.1
    lambda_b.max = max(max(S-diag(p)),-min(S-diag(p)))
    lambda_b.min = lambda_b.min.ratio*lambda_b.max
    lambda_b = exp(seq(log(lambda_b.max), log(lambda_b.min), length = nlambda_b))
  }
  if(is.null(c))
  {
    if(is.null(nc))
      nc = 8
    if(is.null(c.min.ratio))
      c.min.ratio = 0.1
    if(is.null(c.max))
      c.max = 10
    c.min = c.min.ratio*c.max
    c = exp(seq(log(c.max), log(c.min), length = nc))
  }
  
  hpars <- vector(mode = "list", length = 4)
  
  if(is.null(getc)){
    alpha=1/(1+c)
    hpars_table<-expand.grid(alpha, lambda_w, lambda_b)
    key<-hpars_table[, 1]*(hpars_table[, 2]+hpars_table[, 3])
    hpars_table<-cbind(hpars_table, key)
    hpars_table<-hpars_table[order(hpars_table$key, decreasing = T),]
  }
  else{
    hpars_table<-expand.grid(c, lambda_w, lambda_b)
    key<-hpars_table[, 1]*(hpars_table[, 2]+hpars_table[, 3])
    hpars_table<-cbind(hpars_table, key)
    hpars_table<-hpars_table[order(hpars_table$key, decreasing = F),]
    
  }
  
  hpars_table<-hpars_table[,c(1,2,3)]
  hpars[[1]]<-as.matrix(hpars_table)
  if (is.null(getc)) {
    colnames(hpars[[1]])<-c("alpha", "lambda_w", "lambda_b")
  }
  else {
    colnames(hpars[[1]])<-c("c", "lambda_w", "lambda_b")
  }
  rownames(hpars[[1]]) <- seq(1:nrow(hpars[[1]]))
  hpars[[2]]<-lambda_w
  hpars[[3]]<-lambda_b
  hpars[[4]]<-c
  
  hpars
}
