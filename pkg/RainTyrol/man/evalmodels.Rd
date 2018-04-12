\name{evalmodels}
\alias{evalmodels}

\title{Fitting and Comparing Zero-Censored Gaussian Models on Precipitation Data}

\description{
  The function \code{evalmodels} fits distributional trees (\code{\link[disttree]{disttree}}),
  distributional forests (\code{\link[disttree]{distforest}}),
  a prespecified GAMLSS (\code{\link[gamlss]{gamlss}}),
  a boosted GAMLSS (\code{\link[gamboostLSS]{gamboostLSS}}),
  and an EMOS model (\code{\link[crch]{crch}}) to precipitation data.
  The results are compared based on CRPS, log-likelihood and RMSE.
}

\usage{
evalmodels(station, train, test, 
           ntree = 100,
           tree_minsplit = 50, tree_minbucket = 20, tree_mincrit = 0.95,
           forest_minsplit = 50, forest_minbucket = 20, forest_mincrit = 0,
           forest_mtry = 27,
           gamboost_cvr = FALSE)
}

\arguments{
  \item{station}{\code{character}, name of the selected observation station.}
  \item{train}{\code{numeric}, (vector of) years the models should be trained on (available: 1985--2012)}
  \item{test}{\code{numeric}, (vector of) years the models should be tested on (available: 1985--2012)}
  \item{ntree}{\code{numeric}, number of trees in the distributional forest.}
  \item{tree_minsplit}{\code{numeric}, the minimum sum of weights in a node in order to be considered for splitting in the distributional tree.}
  \item{tree_mincrit}{\code{numeric}, the value of the test statistic or 1 - p-value that must be exceeded in order to implement a split in the distributional tree.}
  \item{tree_minbucket}{\code{numeric}, the minimum sum of weights in a terminal node in the distributional tree.}
  \item{forest_minsplit}{\code{numeric}, the minimum sum of weights in a node in order to be considered for splitting in the distributional forest.}
  \item{forest_minbucket}{\code{numeric}, the minimum sum of weights in a terminal node in the distributional forest.}
  \item{forest_mincrit}{\code{numeric}, the value of the test statistic or 1 - p-value that must be exceeded in order to implement a split in the distributional forest.}
  \item{forest_mtry}{\code{numeric}, number of input variables randomly sampled as candidates at each node for random forest like algorithms. The default mtry = Inf means that no random selection takes place.}
  \item{gamboost_cvr}{\code{logical}, Should \code{\link[gamboostLSS]{cvrisk}} be applied to find the optimal value for 'mstop'.}
}

\value{
  \code{evalmodels} returns a list with the following components:
  \item{CRPS}{CRPS (continuos ranked probability score) of all methods, average over testing data.}
  \item{LS}{Logarithmic score (= log-likelihood) of all methods, average over testing data.}
  \item{RMSE}{Root mean squared error of all methods, average over testing data.}
}

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

evalmodels(station = "Axams", train = 1985:2008, test = 2009:2012)

}
}

