\name{dist_vonmises}
\alias{dist_vonmises}
\title{Von Mises Family 'Dist-List' for \code{\link[disttree]{disttree}}.}
\usage{
  dist_vonmises(useC = FALSE, ncores = 1)
}
\arguments{
  \item{useC}{logical; if TRUE C routines are used.}
  \item{ncores}{Number of cores for parallelization with openMP (No big improvements in terms of running time).}
}
\description{
  Exported Von Mises Family for implementation in \code{\link[disttree]{disttree}}.
}
