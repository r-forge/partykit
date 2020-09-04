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
\item{i}{extree data slot index. Can be a numeric vector or \code{"xy"} (see also \code{\link{extree_data}}).
}
  \item{type}{type of variable version to return. 
  \code{"original"} returns the original data;
  \code{"index"} returns the \code{\link[inum]{inum}} object; 
  \code{"scores"} returns the named list of numeric scores (see argument \code{scores} in \code{\link{extree_data}}); 
  \code{"missings"} returns the info on which entries of \code{i} are missing. 
}
}
\details{

}
\value{
The selected information. See argument \code{type}.
}
\seealso{\code{\link{extree_data}}}
\examples{

  data("iris")

  ed <- extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length, 
                    data = iris, nmax = c("yx" = 25, "z" = 10), yx = "matrix")


  ### binned Petal.Width
  extree_variable(ed, 4, type = "index")
  
  ## y and x variables
  head(extree_variable(ed, i = "yx"))
  extree_variable(ed, i = "yx", type = "missings")
  
  ### response
  extree_variable(ed, i = 1)
  
  ### split variables
  extree_variable(ed, i = 4)
  extree_variable(ed, i = 5, type = "index")
  extree_variable(ed, i = 5)[1:10]
  extree_variable(ed, i = c(2, 4))

  
  
  
  
  data("HuntingSpiders", package = "partykit")

  ed_hs <- extree_data(arct.lute + pard.lugu + zora.spin + pard.nigr + 
      pard.pull + aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce + 
      alop.fabr + arct.peri ~ herbs + reft + moss + sand + twigs + water,
      data = HuntingSpiders, yx = "matrix", nmax = c(yx = Inf, z = Inf))
  
  ### the model.frame
  mf <- model.frame(ed)
  all.equal(mf, HuntingSpiders[, names(mf)])
  
  ## Original variable
  extree_variable(ed_hs, 1, type = "original")
  HuntingSpiders$arct.lute
  
  # outcome variables
  head(extree_variable(ed_hs, i = "yx"))
  
  # Several variables
  extree_variable(ed_hs, 10:15)
  
  ## Index 
  # nmax = c(yx = Inf, z = Inf) and thus NULL
  # 1 = arct.lute
  extree_variable(ed_hs, 1, type = "index")
  
  # for z index is created even if nmax = c(yx = Inf, z = Inf)
  # 13 = herbs
  extree_variable(ed_hs, 13, type = "index")
  
  ## No missings
  extree_variable(ed_hs, 1, type = "missings")
  
  
  
  if (require("TH.data") && require("survival")) {
    library("survival")
    data("GBSG2", package = "TH.data")
    
    ed_bc <- extree_data(Surv(time, cens) ~ tgrade, data = GBSG2)
    
    ## Survival outcome
    extree_variable(ed_bc, i = "yx")
    
    ## Ordered variable tgrade (tumor grade)
    extree_variable(ed_bc, i = 2)
    extree_variable(ed_bc, i = 2, type = "scores")
  } 
}

\keyword{tree}
