\name{evalmodels}
\alias{evalmodels}

\title{Fitting and comparing zero-censored Gaussian models on precipitaion data}

\description{
  The function \code{evalmodels} fits distributional trees (\code{disttree}),
  distributional forests (\code{distforest}),
  a prespecified GAMLSS (\code{gamlss}),
  a boosted GAMLSS (\code{gamboostLSS}),
  and an EMOS model (\code{crch}) to precipitation data.
  The results are compared based on CRPS, loglikelihood and RMSE.
}

\usage{
evalmodels(station, train, test, 
           ntree = 100,
           tree_minsplit = 50, tree_minbucket = 20, tree_mincrit = 0.95,
           forest_minsplit = 50, forest_minbucket = 20, forest_mincrit = 0,
           forest_mtry = 27,
           type.tree = "ctree",
           gamboost_cvr = FALSE)
}

\arguments{
  \item{station}{\code{character}, name of the selected observation station.}
  \item{train}{\code{numeric}, (vector of) years the models should be trained on (from 1985--2012)}
  \item{test}{\code{numeric}, (vector of) years the models should be tested on (from 1985--2012)}
  \item{ntree}{\code{numeric}, number of trees in the distributional forest.}
  \item{tree_minsplit}{\code{numeric}, the minimum sum of weights in a node in order to be considered for splitting in the distributional tree.}
  \item{tree_mincrit}{\code{numeric}, the value of the test statistic or 1 - p-value that must be exceeded in order to implement a split in the distributional tree.}
  \item{tree_minbucket}{\code{numeric}, the minimum sum of weights in a terminal node in the distributional tree.}
  \item{forest_minsplit}{\code{numeric}, the minimum sum of weights in a node in order to be considered for splitting in the distributional forest.}
  \item{forest_minbucket}{\code{numeric}, the minimum sum of weights in a terminal node in the distributional forest.}
  \item{forest_mincrit}{\code{numeric}, the value of the test statistic or 1 - p-value that must be exceeded in order to implement a split in the distributional forest.}
  \item{forest_mtry}{\code{numeric}, number of input variables randomly sampled as candidates at each node for random forest like algorithms. The default mtry = Inf means that no random selection takes place.}
  \item{type.tree}{\code{character}, tree algorithm from \pkg{partykit} to be used in \code{disttree}. Can be 'mob' or 'ctree'.}
  \item{gamboost_cvr}{\code{logical}, Should \code{\link[gamboostLSS]{cvrisk}} be applied to find the optimal value for 'mstop'.}
}

\value{
  \code{evalmodels} returns a list with the following components:
  \item{CRPS}{CRPS of all methods, average over testing data}
  \item{LS}{Logarithmic score of all methods, average over testing data}
  \item{RMSE}{Rooted mean squared error of all methods, average over testing data}
}

%\references{}

\examples{
if(require("crch") &
   require("disttree") &
   require("gamlss") &
   require("gamlss.dist") &
   require("gamlss.cens") &
   require("gamboostLSS") &
   require("mboost") &
   require("partykit") &
   require("scoringRules") &
   require("survival")
) {

evalmodels(station = "Axams",
           train = c(1985:2008),
           test = c(2009:2012))

}
}

