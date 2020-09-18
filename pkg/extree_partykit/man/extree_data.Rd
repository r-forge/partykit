\name{extree_data}
\alias{extree_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Data Preprocessing for Extensible Trees.
}
\description{
  A routine for preprocessing data before an extensible tree can be grown by
  \code{extree} or \code{extree_fit}.
}
\usage{
extree_data(formula, data, subset, na.action = na.pass, weights, offset, 
            cluster, strata, scores = NULL, yx = c("none", "matrix"),
            ytype = c("vector", "data.frame", "matrix"), 
            nmax = c(yx = Inf, z = Inf), ...)
}
\arguments{
  \item{formula}{a (multi-part) formula describing the model (\code{y1 + y2 +
		  ... ~ z1 + z2 + ...} or \code{y1 + y2 + ... ~ x1 + x2 + ... |
		  z1 + z2 + ...} or a named list with entries "y" and "z", and
	  optionally "x", "weights", "offset", "cluster", "strata".
}
  \item{data}{an optional data.frame containing the variables in the model.
}
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.
}
  \item{na.action}{a function which indicates what should happen when the
data contain missing values.
}
  \item{weights}{an optional vector of weights.
}
  \item{offset}{an optional offset vector.
}
  \item{cluster}{an optional factor describing clusters. The interpretation
depends on the specific tree algorithm.
}
  \item{strata}{an optional factor describing strata. The interpretation
depends on the specific tree algorithm.
}
  \item{scores}{an optional named list of numeric scores to be assigned to
ordered factors in the \code{z} part of the formula.
}
  \item{yx}{a character indicating if design matrices shall be computed. 
}
  \item{ytype}{a character indicating how response variables shall be
stored.
}
  \item{nmax}{a numeric vector of length two with the maximal number of
bins in the response and \code{x}-part (first element) and the \code{z}
part. Use \code{Inf} to switch-off binning.
}
  \item{\dots}{additional arguments.
}
}
\details{
  This internal functionality will be the basis of implementations of other
  tree algorithms in future versions. Currently, only \code{ctree} relies on
  this function.
  
  We use the following notation:
  \itemize{
      \item \code{y}: response variable(s)
      \item \code{x}: model covariate(s) (e.g. in \code{\link{mob}}, \code{\link{glmtree}}, etc.)
      \item \code{z}: split variable(s)
  }
  
}
\value{An object of class \code{extree_data}. %TODO add mor info here
}
\seealso{\code{\link{extree_variable}}}
\examples{

  
  if(require("mlbench")) {
    
    ## Pima Indians diabetes data
    data("PimaIndiansDiabetes", package = "mlbench")
    
    # with formula
    ed1 <- extree_data(diabetes ~ glucose |  mass + pedigree + age,
        data = PimaIndiansDiabetes)
    
    # without formula (useful for big data sets)
    data_list <- list(
      y = c("diabetes"),
      x = c("glucose"),
      z = c("mass", "pedigree", "age")
    )
    ed2 <- extree_data(formula = data_list, data = PimaIndiansDiabetes)
    
    ## ed1 and ed2 differ, because they are created in different ways.
    ## They still result in the same tree.
    str(ed1$data, 2)
    str(ed2$data, 2)
    
    dim(ed1$data)
    dim(ed2$data)
    
    all.equal(ed1, ed2)
  }

  
  data("iris")

  ed <- extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length, 
                    data = iris, nmax = c("yx" = 25, "z" = 10), yx = "matrix")

  ### the model.frame
  mf <- model.frame(ed)
  all.equal(mf, iris[, names(mf)])

  ### binned y ~ x part
  model.frame(ed, yxonly = TRUE)

  ### binned Petal.Width
  extree_variable(ed, 4, type = "index")

  ### response
  ed$yx$y

  ### model matrix 
  ed$yx$x
  
  
  # Working with missings 
  ## generate example data
  set.seed(2020)
  n <- 15
  X <- rnorm(n = n)
  Z <- rnorm(n = n)
  Y <- rnorm(n = n, mean = x + I(z > 0))
  Zmiss <- Z
  Zmiss[3:5] <- NA
  O <- sample(1:2, size = n, replace = TRUE)
  Omiss <- O
  Omiss[1:2] <- NA
  dat <- data.frame(Y, X, Z, Zmiss, O, Omiss)
    
  exd <- extree_data(formula = list(y = c("Y"), x = c("X"), z = c("Zmiss"), 
      offset = c("Omiss")), data = dat)
    
  ## show all missings that are relevant (not computed for unused variables)
  exd$missings
    
  ## missings in variable Zmiss
  extree_variable(x = exd, i = "Zmiss", type = "missings")
    
  ## missings in yx (this includes the offset, of course)
  extree_variable(x = exd, i = "yx", type = "missings")

}
\keyword{tree}
