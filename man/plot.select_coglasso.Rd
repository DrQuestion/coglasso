% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.select_coglasso}
\alias{plot.select_coglasso}
\alias{plot.coglasso}
\title{Plot selected \code{coglasso} networks}
\usage{
\method{plot}{select_coglasso}(
  x,
  index_c = NULL,
  index_lw = NULL,
  index_lb = NULL,
  node_labels = TRUE,
  hide_isolated = TRUE,
  ...
)

\method{plot}{coglasso}(
  x,
  index_c,
  index_lw,
  index_lb,
  node_labels = TRUE,
  hide_isolated = TRUE,
  ...
)
}
\arguments{
\item{x}{The object of \code{S3} class \code{select_coglasso}.}

\item{index_c}{The index of the \eqn{c} value different from the one
selected by model selection. To set only if the desired network is not the
selected one.}

\item{index_lw}{The index of the \eqn{\lambda_w} value of the chosen
non-optimal network. To set only if the desired network is not the
selected one.}

\item{index_lb}{The index of the \eqn{\lambda_b} value of the chosen
non-optimal network. To set only if the desired network is not the
selected one.}

\item{node_labels}{Show node names in the network. Defaults to TRUE.}

\item{hide_isolated}{Hide nodes that are not connected to any other node.
Defaults to TRUE.}

\item{...}{System required, not used here.}
}
\value{
Returns NULL, invisibly.
}
\description{
\code{plot.select_coglasso()} creates an annotated plot of a \code{coglasso} selected
network from an object of \code{S3} class \code{select_coglasso}. Variables from
different data sets will have different color coding. To plot the network,
it's enough to use \code{plot()} call on the \code{select_coglasso} object.

\code{plot.coglasso()} has the same functioning as \code{select_coglasso.plot()}, but
from an object of \code{S3} class \code{coglasso}. In this case, it is compulsory to
specify \code{index_c}, \code{index_lw}, and \code{index_lb.}
}
\details{
If the input is a \code{coglasso} object, it is necessary to specify all the
indexes to extract a selected network. \cr
If the input is a \code{select_coglasso} object, it extracts by default the
selected network. If the selection method was "ebic", and you want to extract
a different network than the selected one, specify all indexes.
Otherwise, if the objective is to extract the optimal network for a specific
\eqn{c} value different than the selected one, set \code{index_c} to your chosen
one. Also here it is possible to extract a specific non-optimal network by
setting all the indexes to the chosen ones.
}
\examples{

\donttest{
sel_cg <- bs(multi_omics_sd_small, p = c(14, 5), nlambda_w = 15, nlambda_b = 15,
             nc = 3, lambda_w_min_ratio = 0.6, verbose = FALSE)
plot(sel_cg)
}

}
\seealso{
\code{\link[=get_network]{get_network()}} to understand what it means to select a specific
network with \code{index_c}, \code{index_lw}, and \code{index_lb.}
}
