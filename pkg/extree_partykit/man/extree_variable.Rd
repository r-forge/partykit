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
extree_variable(x, index, variable, type = c("original", "inum", "scores", "missings"))
}
\arguments{
  \item{x}{an object of class \code{\link{extree_data}}.
}
  \item{index}{extree data slot index. Can be a numeric vector or a character vector of variable names.
  If \code{NULL}, \code{variable} must be specified. 
}
  \item{variable}{extree data slot name. Can be \code{"yx"} or \code{"y"} (see also \code{\link{extree_data}}).
  If \code{NULL}, \code{index} must be specified. 
}
  \item{type}{type of variable version to return. 
  \code{"original"} returns the original data;
  \code{"inum"} returns the \code{\link[inum]{inum}} object; 
  \code{"scores"} returns the named list of numeric scores (see argument \code{scores} in \code{\link{extree_data}}); 
  \code{"missings"} returns the info on which entries of \code{index} are missing. 
}
}
\details{

}
\value{
The selected information. See argument \code{type}.
}
\seealso{\code{\link{extree_data}}}
\examples{

data("iris", package = "datasets")
ed <- extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length, 
data = iris, nmax = c("yx" = 25, "z" = 10), yx = "matrix")

### binned Petal.Width
extree_variable(ed, index = 4, type = "inum")

## y and x variables
head(extree_variable(ed, variable = "yx"))
extree_variable(ed, variable = "yx", type = "missings")

### response
extree_variable(ed, index = 1)
extree_variable(ed, index = "Species")

### split variables
extree_variable(ed, index = 4)
extree_variable(ed, index = "Petal.Width")
extree_variable(ed, index = 5, type = "inum")
extree_variable(ed, index = 5)[1:10]

### several variables at once
## type = "original"
extree_variable(ed, index = c(2, 4))
extree_variable(ed, index = c("Sepal.Width", "Petal.Width"))
## type = "inum" (only for binned variables, otherwise NULL)
extree_variable(ed, index = c(2, 4), type = "inum")
extree_variable(ed, index = c("Sepal.Width", "Petal.Width"), type = "inum")


data("HuntingSpiders", package = "partykit")

ed_hs <- extree_data(arct.lute + pard.lugu + zora.spin + pard.nigr + 
    pard.pull + aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce + 
    alop.fabr + arct.peri ~ herbs + reft + moss + sand + twigs + water,
  data = HuntingSpiders, yx = "matrix", nmax = c(yx = Inf, z = Inf))

### the model.frame
mf <- model.frame(ed_hs)
## mf differs from HuntingSpiders because of its additional attributes 
## derived from its formula, although the data.frames are the same
all.equal(mf, HuntingSpiders[, names(mf)])

## Original variable
extree_variable(ed_hs, 1, type = "original")
HuntingSpiders$arct.lute

# outcome variables
head(extree_variable(ed_hs, variable = "yx"))

# Several variables
extree_variable(ed_hs, 10:15)

## Index 
# nmax = c(yx = Inf, z = Inf) and thus NULL
# 1 = arct.lute
extree_variable(ed_hs, 1, type = "inum")

# for z index is created even if nmax = c(yx = Inf, z = Inf)
# 13 = herbs
extree_variable(ed_hs, 13, type = "inum")

## No missings
extree_variable(ed_hs, 1, type = "missings")



if (require("TH.data") && require("survival")) {
  library("survival")
  data("GBSG2", package = "TH.data")
  
  ed_bc <- extree_data(Surv(time, cens) ~ tgrade, data = GBSG2)
  
  ## Survival outcome
  extree_variable(ed_bc, variable = "yx")
  
  ## Ordered variable tgrade (tumor grade)
  extree_variable(ed_bc, index = 2)
  extree_variable(ed_bc, index = 2, type = "scores")
} 
}

\keyword{tree}
