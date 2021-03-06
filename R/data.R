#' Comparative dataset used to explore the macroevolution of echinoid body size.
#'
#' A list containing a data.set of body sizes, a vector of measurement errors, a
#' time-calibrated phylogeny, and the results of running forward and backward
#' phases of SURFACE on this dataset.
#'
#' @format A list with all the data to replicate the macroevolutionary analysis
#'   in Mongiardino Koch & Thompson (2020):
#' \describe{
#'   \item{size}{A data.frame including body sizes for all terminals and
#'   terminal names as row names}
#'   \item{error}{A named vector with measurement errors, in the same order as \code{$size}}
#'   \item{tree}{A time-calibrated phylogeny of all taxa with size data}
#'   \item{fwd_surface}{The result of running \code{surfaceForward} on this dataset}
#'   \item{bwd_surface}{The result of running \code{surfaceBackward} on this dataset}
#' }
#' @references Mongiardino Koch N., Thompson J.R. A Total-Evidence Dated
#'   Phylogeny of Echinoids and the Evolution of Body Size across Adaptive
#'   Landscape. bioRxiv 2020.02.13.947796, doi.org/10.1101/2020.02.13.947796
"echinoid_data"
