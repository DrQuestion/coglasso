#' Plot selected `coglasso` networks
#' 
#' @description
#' `plot.select_coglasso()` creates an annotated plot of a `coglasso` selected 
#' network from an object of `S3` class `select_coglasso`. Variables from 
#' different data sets will have different color coding. To plot the network, 
#' it's enough to use `plot()` call on the `select_coglasso` object.
#' 
#' `plot.coglasso()` has the same functioning as `select_coglasso.plot()`, but 
#' from an object of `S3` class `coglasso`. In this case, it is compulsory to 
#' specify `index_c`, `index_lw`, and `index_lb.`
#'
#' @param x The object of `S3` class `select_coglasso`.
#' @inherit get_network
#' @param node_labels Show node names in the network. Defaults to TRUE. 
#' @param hide_isolated Hide nodes that are not connected to any other node.
#'   Defaults to TRUE.
#' @param ... System required, not used here.
#'
#' @seealso [get_network()] to understand what it means to select a specific 
#' network with `index_c`, `index_lw`, and `index_lb.`
#' 
#' @return Returns NULL, invisibly.
#' @export
#'
#' @examples
#' 
#' \donttest{
#' sel_cg <- bs(multi_omics_sd_small, p = c(14, 5), nlambda_w = 15, nlambda_b = 15,
#'              nc = 3, lambda_w_min_ratio = 0.6, verbose = FALSE)
#' plot(sel_cg)
#' }
#' 
plot.select_coglasso <- function(x, index_c=NULL, index_lw=NULL, index_lb=NULL, node_labels = TRUE, hide_isolated = TRUE, ...) {
  sel_network <- get_network(x, index_c=index_c, index_lw=index_lw, index_lb=index_lb)
  sel_pcor <- get_pcor(x, index_c=index_c, index_lw=index_lw, index_lb=index_lb)
  
  if (!node_labels) {
    igraph::V(sel_network)$label <- NA
  }
  
  p <- x$p
  p_tot <- ncol(x$data)
  D <- x$D
  fillsframes <- get_fillsframes(D)
  igraph::V(sel_network)$color <- rep(0, p_tot)
  igraph::V(sel_network)$frame.color <- rep(0, p_tot)
  igraph::V(sel_network)$label.color <- rep(0, p_tot)
  i <- 1
  for (j in 1:p_tot){
    if (j > sum(p[1:i])) {
      i <- i + 1
    }
    igraph::V(sel_network)$color[j] <- fillsframes[[1]][i]
    igraph::V(sel_network)$frame.color[j] <- fillsframes[[2]][i]
    igraph::V(sel_network)$label.color[j] <- fillsframes[[2]][i]
  }
  
  ws <- abs(sel_pcor)
  upper <- c(ws[upper.tri(ws)])
  lower <- c(t(ws)[upper.tri(ws)])
  ws <- sapply(seq_along(upper), function(i) max(upper[i], lower[i]))
  ws <- ws[ws != 0]
  
  lo <- igraph::layout_with_fr(sel_network, weights = ws)
  if (hide_isolated) {
    disconnected <- which(igraph::degree(sel_network) == 0)
    sel_network <- igraph::delete_vertices(sel_network, disconnected)
    if (length(disconnected)) {
      lo <- lo[-disconnected, ]
    }
  }
  
  igraph::V(sel_network)$frame.width <- 5.5/log(igraph::gorder(sel_network)+1)
  igraph::V(sel_network)$size <- 40/log(igraph::gorder(sel_network)+1)
  igraph::E(sel_network)$width <- 4/log(igraph::gsize(sel_network)+1)
  
  plot(sel_network, layout = lo)
  return(invisible(NULL))
}

#' @rdname plot.select_coglasso
#' @export
plot.coglasso <- function(x, index_c, index_lw, index_lb, node_labels = TRUE, hide_isolated = TRUE, ...) {
  sel_network <- get_network(x, index_c=index_c, index_lw=index_lw, index_lb=index_lb)
  sel_pcor <- get_pcor(x, index_c=index_c, index_lw=index_lw, index_lb=index_lb)
  
  if (!node_labels) {
    igraph::V(sel_network)$label <- NA
  }
  
  p <- x$p
  p_tot <- ncol(x$data)
  D <- x$D
  fillsframes <- get_fillsframes(D)
  igraph::V(sel_network)$color <- rep(0, p_tot)
  igraph::V(sel_network)$frame.color <- rep(0, p_tot)
  igraph::V(sel_network)$label.color <- rep(0, p_tot)
  i <- 1
  for (j in 1:p_tot){
    if (j > sum(p[1:i])) {
      i <- i + 1
    }
    igraph::V(sel_network)$color[j] <- fillsframes[[1]][i]
    igraph::V(sel_network)$frame.color[j] <- fillsframes[[2]][i]
    igraph::V(sel_network)$label.color[j] <- fillsframes[[2]][i]
  }
  
  ws <- abs(sel_pcor)
  upper <- c(ws[upper.tri(ws)])
  lower <- c(t(ws)[upper.tri(ws)])
  max_ws <- sapply(seq_along(upper), function(i) max(upper[i], lower[i]))
  ws[upper.tri(ws, diag = FALSE)] <- max_ws
  ws <- t(ws)
  ws[upper.tri(ws, diag = FALSE)] <- max_ws
  igraph::E(sel_network)$weight <- ws[igraph::as_edgelist(sel_network)]
  
  lo <- igraph::layout_with_fr(sel_network)
  if (hide_isolated) {
    disconnected <- which(igraph::degree(sel_network) == 0)
    sel_network <- igraph::delete_vertices(sel_network, disconnected)
    if (length(disconnected)) {
      lo <- lo[-disconnected, ]
    }
  }
  
  igraph::V(sel_network)$frame.width <- 5.5/log(igraph::gorder(sel_network)+1)
  igraph::V(sel_network)$size <- 40/log(igraph::gorder(sel_network)+1)
  igraph::E(sel_network)$width <- 4/log(igraph::gsize(sel_network)+1)
  
  plot(sel_network, layout = lo)
  return(invisible(NULL))
}

#' Extract colors for nodes and frames to color code variables from different data sets
#'
#' @param num_datasets The number of data sets. Up to six, the returned colors 
#'   are part of a curated palette.
#'
#' @return `get_fillsframes()` returns a list of two paired vectors, with as 
#'   many elements as `num_datasets`, one with node colors and the other one 
#'   with the relative darker node frame colors.
#'   
#' @noRd
#' 
get_fillsframes <- function(num_datasets){
  fills <- c("#00ccff", "#ff9999", "#F2F79E", "#7AD988", "#AFA2FF", "#A63A50")
  frames <- c("#002060", "#800000", "#66684d", "#0c3b12", "#48417c", "#5e2d38")
  if (num_datasets <= length(fills)) {
    return(list(fills[1:num_datasets], frames[1:num_datasets]))
  }
  else{
    add_fills <- grDevices::rainbow(num_datasets-length(fills))
    add_frames <- sapply(add_fills, get_frame, USE.NAMES = FALSE)
    return(list(c(fills, add_fills), c(frames, add_frames)))
  }
}

#' Get a frame color relative to a given node color
#'
#' @param color The node color.
#'
#' @return Returns the frame color relative to the input node color, a tint.
#'   
#' @noRd
#' 
get_frame <- function(color){
  grDevices::colorRampPalette(colors = c(color, "black"))(11)[8]
}
