library("zoo")

# load and prepare data
d <- readRDS("~/svn/partykit/pkg/disttree/inst/circforest/circforest_data/circforest_prepared_data_ibk_lag3.rds")
dim(d)

#nasums <- rep.int(0,NCOL(dz))
#for(i in 1:NCOL(dz)) nasums[i] <- sum(is.na(dz[,i]))
#dzid <- which(!is.na(dz$dd.response))
#dz <- dz[dzid,]

d <- fortify.zoo(d)
head(d[,1:5])
dim(d)
#dall <- d
d <- d[complete.cases(d),]
dim(d)
names(d)
head(d[,1:5])

# transform response.dd from 0-360 degree to 0-2pi
summary(d$dd.response)
d$dd.response <- d$dd.response / 360 * 2*pi

## save reduced data set
save(d, file = "~/svn/partykit/pkg/disttree/inst/circforest/circforest_data/d_prep.rda")

########################################3
# tests 

load("~/svn/partykit/pkg/disttree/inst/circforest/circforest_data/d_prep.rda")

## FIXME: dist_vonmises operates on [-pi, pi]
## transform dd.response to [-pi,pi) but without transforming 0, 
## hence, only transform [pi, 2*pi) to [-pi,0)
d$dd.response[d$dd.response>pi] <- d$dd.response[d$dd.response>pi] - 2*pi
summary(d$dd.response)
library("disttree")
library("circmax")
source("~/svn/uibk-rprog-2018/pkg/circmax/R/vonmises.R")

f <- as.formula(paste("dd.response ~ ", paste(names(d)[-c(1:8)], collapse= "+")))

## fit to full data set
dt <- distextree(formula = f,
                 data = d, 
                 family = dist_vonmises,
                 control = distextree_control(maxdepth = 4)) 
plot(dt)

## FIXME: minbucket and minsplit have to be set -> set default values
# (otherwise error in .extree_node where minbucket/minsplit is compared to nr of current values)
df <- distexforest(formula = f,
                   data = d, 
                   family = dist_vonmises,
                   ntree = 100,
                   perturb = list(replace = FALSE, fraction = 0.632),
                   control = distextree_control(minbucket = 1000, minsplit = 1500))

predict(df, newdata = d[c(100),])
predict(df, newdata = d[c(100),], type = "response")   ## FIXME: caculation of expected value (integrate)
predict(df, newdata = d[c(100),], type = "node")


## FIXME: memory error: cannot allocate vector of size 119.7 Gb
vi <- disttree:::varimp.distexforest(df, nperm = 10)


## evaluate performance with circular crps

# first define learndata and testdata
summary(d$Index)
NROW(d[d$Index > "2017-12-31 00:00:00",])
NROW(d[d$Index > "2016-12-31 00:00:00" & d$Index < "2017-12-31 00:00:00",])
NROW(d[d$Index > "2015-12-31 00:00:00" & d$Index < "2016-12-31 00:00:00",])
NROW(d[d$Index > "2014-12-31 00:00:00" & d$Index < "2015-12-31 00:00:00",])
NROW(d[d$Index > "2013-12-31 00:00:00" & d$Index < "2014-12-31 00:00:00",])

dlearn <- d[d$Index <= "2018-01-01 00:00:00" & d$Index > "2015-01-01 00:00:00",]
dtest <- d[d$Index > "2018-01-01 00:00:00",]

dt <- distextree(formula = f,
                 data = dlearn, 
                 family = dist_vonmises,
                 control = distextree_control(maxdepth = 4)) 
plot(dt)
predpar_dt <- predict(dt, newdata = dtest, type = "parameter")
crps_dt <- circmax:::crps_vonmises(dtest$dd.response, 
                                   mu = predpar_dt$mu, 
                                   kappa = predpar_dt$kappa,
                                   sum = TRUE)


set.seed(7)
df <- distexforest(formula = f,
                   data = dlearn, 
                   family = dist_vonmises,
                   ntree = 15,
                   perturb = list(replace = FALSE, fraction = 0.5),
                   control = distextree_control(minbucket = 1000, minsplit = 4000))

# predict parameter stepwise (due to memory problems)
ntest <- NROW(dtest)
n10 <- round(ntest/10)
predpar_df <- predict(df, newdata = dtest[c(1:n10),], type = "parameter")
for(i in 1:8){
  predpar_df_part <- predict(df, newdata = dtest[c((n10*i + 1) : (n10*(i+1))),], type = "parameter")
  predpar_df <- rbind(predpar_df,predpar_df_part)
}
predpar_df_part <- predict(df, newdata = dtest[c((n10*9 + 1):ntest),], type = "parameter")
predpar_df <- rbind(predpar_df,predpar_df_part)

crps_df <- circmax:::crps_vonmises(dtest$dd.response, 
                                   mu = predpar_df$mu, 
                                   kappa = predpar_df$kappa,
                                   sum = TRUE)

crps <- c(crps_dt, crps_df)
names(crps) <- c("crps_dt", "crps_df")
crps



## climatology
summary(d$Index)

dim(dlearn)
dim(dtest)
nt <- NROW(dtest)
mu <- kappa <- numeric(length = nt)

for(i in 1:nt){
  
  id <- dtest$Index[i]
  
  learnids_hour_each_year <- c(id - 60*60*24*365 + 10*60 * c(-3:3),
                               id - 60*60*24*(365+366) + 10*60 * c(-3:3),
                               id - 60*60*24*(365*2+366) + 10*60 * c(-3:3))
  
  learnids <- learnids_hour_each_year[1] + + 60*60*24 * c(-10:10)
  for(j in 2:length(learnids_hour_each_year)){
    learnids_10days <- learnids_hour_each_year[j] + + 60*60*24 * c(-10:10)
    learnids <- c(learnids, learnids_10days)
  }
  
  mu[i] <- mean(dlearn$dd.response[dlearn$Index %in% learnids]) 
  kappa[i] <- 1/var(dlearn$dd.response[dlearn$Index %in% learnids]) 
                
}

climpar <- data.frame("mu" = mu,
                      "kappa" = kappa)

crps_clim <- circmax:::crps_vonmises(dtest$dd.response, 
                                     mu = climpar$mu, 
                                     kappa = climpar$kappa,
                                     sum = TRUE)

crps <- c(crps_dt, crps_df, crps_clim)
names(crps) <- c("crps_dt", "crps_df", "crps_clim")
crps

head(cbind(climpar$mu, predpar_df$mu, predpar_dt$mu, dtest$dd.response))






if(FALSE){
  
  
  # try without differences and change variables
  diff_ids <- grepl("diff", names(d))
  change_ids <- grepl("_ch", names(d))
  cn <- names(d)[!(diff_ids | change_ids)]
  cn <- cn[-c(1:8)]
  length(cn)
  
  f <- as.formula(paste("dd.response ~ ", paste(cn, collapse= "+")))
  
  dt <- distextree(formula = f,
                   data = d, 
                   family = dist_vonmises,
                   control = distextree_control(minsplit = 10000)) # alpha = 0.05 per default(?)
  plot(dt)
  ## FIXME: does not detect any split unless alpha is set to a higher values than 0.05
  
  dt <- distextree(formula = f,
                   data = d, 
                   family = dist_vonmises,
                   control = distextree_control(minsplit = 10000, mincriterion = 0.2)) # alpha = 0.8
  plot(dt)
  
  
  ## trying different variables
  f <- as.formula(paste("dd.response ~ ", 
                        paste(colnames(d[,c(9:17,19:20,24:26,30:38,43,44,51,
                                            177:195,202:210, 214:234, 241)]), 
                              collapse= "+")))
  
  dt <- distextree(formula = f,
                   data = d, 
                   family = dist_vonmises,
                   control = distextree_control(minsplit = 10000, 
                                                tol = sqrt(.Machine$double.eps)^1))
  
  plot(dt)
  names(d)[c(18,21:23,27:29,39:42,45:50,
             196:201, 211:213,235:240)]
  
}