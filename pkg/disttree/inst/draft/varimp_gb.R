###################
# Varimp for gamboostLSS
###################

nperm <- 50
  
library("gamboostLSS")
library("crch")
library("scoringRules")
library("RainTyrol")

# set function for parallelization
applyfun <- function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = pmax(1, detectCores() - 1))


# model
load("~/svn/partykit/pkg/disttree/inst/draft/Axams_24to4.rda")
gb <- Axams_24to4$gb

# data
data("RainTyrol")
RainData <- RainTyrol[RainTyrol$station == "Axams", ]
rownames(RainData) <- c(1:NROW(RainData))
  
# learning data: 24 years (1985 - 2008, both inlcuded)
# testing data: 4 successive years (2009, 2010, 2011, 2012)
# learndata <- RainData[RainData$year < 2009,]
testdata <- RainData[RainData$year %in% c(2009, 2010, 2011, 2012),]
  
  
## compute variable importance for the gamboostLSS object 'gb'

set.seed(7)

meancrps <- function(permute = NULL, newdata = testdata) {
  if(!is.null(permute)) newdata[[permute]] <- sample(newdata[[permute]])
  m <- predict(gb, newdata = newdata, parameter = "mu", type = "response")
  s <- predict(gb, newdata = newdata, parameter = "sigma", type = "response")
  mean(crps_cnorm(newdata$robs, location = m, scale = s, lower = 0))
}


# apply for all covariates except for dswrf_mean_min and 
# dswrf_sprd_min (columns 30 and 33) as they are always 0

# risk <- sapply(c(5:29, 31, 32, 34: ncol(testdata)), meancrps)

# using only one core
# risk_all <- replicate(nperm, sapply(c(5:29, 31, 32, 34: ncol(testdata)), meancrps))
# risk <- rowMeans(risk_all)

# or parallel
risklist <- applyfun(1:nperm, 
                     function(i){
                       sapply(c(5:29, 31, 32, 34: ncol(testdata)), meancrps)
                     })
risk <- Reduce("+", risklist) / length(risklist)

names(risk) <- names(testdata)[c(5:29, 31, 32, 34: ncol(testdata))]
vimp_crps_gb <- risk - meancrps(newdata = testdata)
vimp_crps_gb <- sort(vimp_crps_gb, decreasing = TRUE)

save(vimp_crps_gb, file = "~/svn/partykit/pkg/disttree/inst/draft/reviews/vimp_crps_gb.rda")

top10_gb <- vimp_crps_gb[1:10]

par(mar = c(2.5,10,1,2))
barplot(rev(top10_gb), 
        horiz = TRUE, las = 1, axes = FALSE,
        # xlab = "Variable importance: mean decrease in CRPS",
        font.axis = 3, #list(family="HersheySerif", face=3),
        # xlim = c(0,1.1),
        xlim = c(0,0.4),
        names.arg = gsub("pow", "", names(rev(top10_gb))))
axis(1, at = seq(0,1.6,0.1), las = 1, mgp=c(0,1,0))