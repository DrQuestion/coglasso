#' Extract a `coglasso` network
#' 
#' `get_network()` extracts the selected network from a `select_coglasso` object,
#' or a different specific one from either a `select_coglasso` or a `coglasso` 
#' object when specifying the optional parameters. 
#' 
#' If the input is a `coglasso` object, it is necessary to specify all the
#' indexes to extract a selected network. \cr
#' If the input is a `select_coglasso` object, it extracts by default the 
#' selected network. If the selection method was "ebic", and you want to extract 
#' a different network than the selected one, specify all indexes. 
#' Otherwise, if the objective is to extract the optimal network for a specific 
#' \eqn{c} value different than the selected one, set `index_c` to your chosen 
#' one. Also here it is possible to extract a specific non-optimal network by 
#' setting all the indexes to the chosen ones.
#' 
#' @encoding UTF-8
#' @param sel_cg_obj The object of `S3` class `select_coglasso` or of `S3` 
#'   class `coglasso`.
#' @param index_c The index of the \eqn{c} value different from the one 
#'   selected by model selection. To set only if the desired network is not the 
#'   selected one.
#' @param index_lw The index of the \eqn{\lambda_w} value of the chosen 
#'   non-optimal network. To set only if the desired network is not the 
#'   selected one.
#' @param index_lb The index of the \eqn{\lambda_b} value of the chosen 
#'   non-optimal network. To set only if the desired network is not the 
#'   selected one.
#'
#' @return `get_network()` returns the selected network, in the form of an 
#' object of class `igraph`.
#' @export
#'
#' @examples
#' \donttest{
#' sel_cg <- bs(multi_omics_sd_micro, p = c(4, 2), nlambda_w = 3, nlambda_b = 3,
#'                  nc = 3, verbose = FALSE)
#' sel_net <- get_network(sel_cg)
#' # Could even plot the selected network with plot(sel_net), but then it would
#' # plot an unnotated network, better to directly plot(sel_cg).
#' print(sel_net)
#' }
#' 
get_network <- function(sel_cg_obj, index_c=NULL, index_lw=NULL, index_lb=NULL){
  # Dealing with class "coglasso"
  if (inherits(sel_cg_obj, "coglasso") & (is.null(index_lw) | is.null(index_lb) | is.null(index_c))) {
    stop('If input object is of class "coglasso" please select which network to extract. See help(get_network).')
  }
  else if (inherits(sel_cg_obj, "coglasso")) {
    c <- sel_cg_obj$c[index_c]
    D <- sel_cg_obj$D
    alpha <- 1 / (c * (D - 1) + 1)
    index_network <- which(sel_cg_obj$hpars[, 1] == alpha &
                             sel_cg_obj$hpars[, 2] == sel_cg_obj$lambda_w[index_lw] &
                             sel_cg_obj$hpars[, 3] == sel_cg_obj$lambda_b[index_lb] &
                             sel_cg_obj$hpars[, 4] == c)
    network <- igraph::graph.adjacency(sel_cg_obj$path[[index_network]], mode = "undirected")
    if (colnames(sel_cg_obj$data)) {
      igraph::V(network)$label <- colnames(sel_cg_obj$data)
    }
    return(network)
  }
  # Dealing with selection method "ebic"
  if (sel_cg_obj$method=="ebic" & is.null(index_lw) & is.null(index_lb) & is.null(index_c)){
    network <- igraph::graph.adjacency(sel_cg_obj$sel_adj, mode = "undirected")
    if (colnames(sel_cg_obj$data)) {
      igraph::V(network)$label <- colnames(sel_cg_obj$data)
    }
    return(network)
  }
  else if (sel_cg_obj$method=="ebic" & (is.null(index_lw) | is.null(index_lb) | is.null(index_c))) {
    stop('If an index_lw, an index_lb, or an index_c is specified, the others must be specified too.')
  }
  else if (sel_cg_obj$method=="ebic") {
    c <- sel_cg_obj$c[index_c]
    D <- sel_cg_obj$D
    alpha <- 1 / (c * (D - 1) + 1)
    index_network <- which(sel_cg_obj$hpars[, 1] == alpha &
                             sel_cg_obj$hpars[, 2] == sel_cg_obj$lambda_w[index_lw] &
                             sel_cg_obj$hpars[, 3] == sel_cg_obj$lambda_b[index_lb] &
                             sel_cg_obj$hpars[, 4] == c)
    network <- igraph::graph.adjacency(sel_cg_obj$path[[index_network]], mode = "undirected")
    if (colnames(sel_cg_obj$data)) {
      igraph::V(network)$label <- colnames(sel_cg_obj$data)
    }
    return(network)
  }
  # Dealing with remaining cases
  if (!is.null(index_c)){
    if (is.null(index_lw) & is.null(index_lb)) {
      index_lw <- sel_cg_obj$opt_index_lw[index_c]
      index_lb <- sel_cg_obj$opt_index_lb[index_c]
    }
    else if (is.null(index_lw) | is.null(index_lb)) {
      stop('If an index_lw or an index_lb is specified, the other must be specified too.')
    }
    c <- sel_cg_obj$c[index_c]
    D <- sel_cg_obj$D
    alpha <- 1 / (c * (D - 1) + 1)
    index_network <- which(sel_cg_obj$hpars[, 1] == alpha &
                             sel_cg_obj$hpars[, 2] == sel_cg_obj$lambda_w[index_lw] &
                             sel_cg_obj$hpars[, 3] == sel_cg_obj$lambda_b[index_lb] &
                             sel_cg_obj$hpars[, 4] == c)
    network <- igraph::graph.adjacency(sel_cg_obj$path[[index_network]], mode = "undirected")
    if (colnames(sel_cg_obj$data)) {
      igraph::V(network)$label <- colnames(sel_cg_obj$data)
    }
    return(network)
  }
  network <- igraph::graph.adjacency(sel_cg_obj$sel_adj, mode = "undirected")
  if (colnames(sel_cg_obj$data)) {
    igraph::V(network)$label <- colnames(sel_cg_obj$data)
  }
  return(network)
}