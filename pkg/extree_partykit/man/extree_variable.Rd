\name{extree_variable}
\alias{extree_variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract variable(s) from \code{extree_data} object
}
\description{
  \code{extree_data} objects contain a lot of useful information. This convenience
  function helps accessing them.
}
\usage{
extree_variable(x, i, type = c("original", "index", "scores", "missings"))
}
\arguments{
  \item{x}{an object of class \code{\link{extree_data}}.
}
\item{i}{extree data slot index. Can be a numeric vector or \code{"xy"}.
}
  \item{type}{type of variable version to return. 
  \code{"original"} returns the original data;
  \code{"index"} returns the \code{\link[inum]{inum}} object; 
  \code{"scores"} returns the named list of numeric scores (see argument \code{scores} in \code{\link{extree_data}}); 
  \code{"missings"} returns the info on which entries of \code{i} are missing; 
}
}
\details{

}
\value{
The selected information.
}
\keyword{tree}
