% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ttt.R
\name{plot_ttt}
\alias{plot_ttt}
\title{Plot trait dynamics through time}
\usage{
plot_ttt(
  ttt,
  interval = NA,
  window_size = NA,
  CI = 80,
  trim = T,
  graphs = "both",
  rates = "diversification",
  k = NA
)
}
\arguments{
\item{ttt}{Data.frame output by \code{\link{transitions_through_time}}.}

\item{interval}{Numeric value that determines the temporal resolution (i.e.,
the step size in Ma at which the number of active regimes is recorded).
Defaults to slightly less than 100 given the time spanned by the phylogeny.}

\item{window_size}{Numeric value that sets the width of the window used to
smoot rate estimates. Defaults to a width that includes approx. 10
intervals (see above).}

\item{CI}{Numeric value that sets the confidence interval (expressed as
percentage). Determines the amount of results that are discarded before
plotting (default = 80).}

\item{trim}{Whether to trim a few values at the begining and end of plot that
contain fewer intervals and can be noisier. Default is \code{TRUE}.}

\item{graphs}{Which graphs to plot. Options include \code{'active_regimes'},
\code{'rate_ttt'}, and \code{'both'}.}

\item{rates}{Which rates to plot. Options include any combination of \code{'birth'},
\code{'death'}, and \code{'diversification'}. Defaults to only the latter.}

\item{k}{The value of k used for gam regression. If not specified this is
automatically determined (see more details in \code{\link[mgcv]{gam}}).}
}
\value{
A plot including different visual summaries of the evolutionary
  dynamics of discrete traits through time.
}
\description{
This function plots evolutionary dynamics of discrete traits as summarized
using \code{\link{transitions_through_time}}. Two types of plots can be
generated: the number of different states active through time (where identical
states that have different origins are not counted as the same one), and the
rates at which these states are originating (birth rate), becoming extinct
(death rate), or accumulating (diversification rate) thorugh time. Depending
on the character being investigated, these plots might or might not be
meaningful.
}
\details{
By default, this is used by \code{\link{transitions_through_time}} to plot
results. However, the object returned by that function can also be used here
with more control on the plotting options. These include the intervals (in Ma)
at which the number of states are recorded, the size of the window used to
smooth rates, the type of plot generated, and the type of rate to plot (see
Arguments).

Trends are depicted using GAM regressions (see \code{\link[mgcv]{gam}}).
Depending on the combination of the size of the smoothing window and the
number of smoothing functions used, nonsensical results can be obtained. Some
tuning might be necessary to correctly depict trends in the data.
}
\references{
Mongiardino Koch N. 2021. Exploring adaptive landscapes across deep time:
  A case study using echinoid body size. Evolution,
  https://doi.org/10.1111/evo.14219.
}
\seealso{
\code{\link{transitions_through_time}}
}
\author{
Nicolás Mongiardino Koch
}
