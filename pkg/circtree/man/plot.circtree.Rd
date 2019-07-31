\name{plot.circtree}
\alias{plot.circtree}
\title{Plotting a Regression Tree with a Circular Response (under development).}
\usage{
  \method{plot}{circtree}(x, terminal_panel = node_circular,
    tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
}
\arguments{
\item{x}{Object of class \code{\link{circtree}}.}

\item{terminal_panel}{Do not change.}

\item{tp_args}{Do not change.}

\item{tnex}{Do not change.}

\item{drop_terminal}{Do not change.}

\item{\dots}{Do not change.}
}
\description{
  This function plots regression trees with a circular response based on \code{\link[partykit]{plot.constparty}}.
}
