\name{extree_data-methods}
\alias{extree_data-methods}
\alias{dim.extree_data}
\alias{dimnames.extree_data}
\alias{print.extree_data}
\alias{str.extree_data}

\title{ Methods for extree_data objects }
\description{

  Methods for computing on \code{extree_data} objects.

}
\usage{
\method{str}{extree_data}(object, max.level = 1, give.attr = FALSE, ...)
\method{print}{extree_data}(x, ...)
\method{dimnames}{extree_data}(x)
\method{dim}{extree_data}(x)
}
\arguments{
  \item{x}{ an object of class \code{\link{extree_data}}.}
  \item{object}{ an object of class \code{\link{extree_data}}.}
  \item{max.level}{maximal level of nesting which is applied for displaying 
    nested structures, e.g., a list containing sub lists. Default NA: Display 
    all nesting levels.}
  \item{give.attr}{logical; if TRUE (default), show attributes as sub structures.}
  \item{\dots}{ additional arguments passed to \code{\link{str}}.}
}
\details{

  \code{str} compactly shows the internal structure of the object and 
  \code{print} returns dimension of the \code{data.frame} and the names of the variables.
  
  \code{dimnames} gets the \code{dimnames} of \code{x$data} and \code{dim} gets 
  the \code{dim} of \code{x$data}.

}
\seealso{\code{\link{extree_data}}, \code{\link{extree_variable}}}
\examples{

    ## airquality data
    airq <- subset(airquality, !is.na(Ozone))
    my_data <- extree_data(Ozone ~ Wind + Temp, 
        data = airq, yx = "matrix")
        
    ## use methods
    print(my_data)
    str(my_data, 2)
    dim(my_data)
    dimnames(my_data)
    
}
