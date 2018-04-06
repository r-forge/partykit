library("disttree")
library("crch")
library("countreg")
library("gamlss")
library("gamlss.cens")
gen.cens(NO, type = "left")
library("RainTyrol")

## HCL palette
pal <- hcl(c(10, 128, 260, 290, 50), 100, 50)
names(pal) <- c("forest", "tree", "gamlss", "gamboostLSS", "EMOS")

pallight <- hcl(c(10, 128, 260, 290, 70), 100, 50, alpha = 0.25)
names(pallight) <- c("forest", "tree", "gamlss", "gamboostLSS", "EMOS")

transpgray <- rgb(0.190,0.190,0.190, alpha = 0.2)

## set function for parallelization
applyfun <- function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = detectCores() - 1)


#####################################################  
#### cross validation at station Axams

nrep_cross <- 10
seed <- 7

res_cross <- applyfun(1:nrep_cross,
                      function(i){
                        
                        set.seed(seed*i)
                        
                        # randomly split data in 7 parts each including 4 years
                        years <- 1985:2012
                        testyears <- list()
                        for(j in 1:7){
                          testyears[[j]] <- sample(years, 4, replace = FALSE)
                          years <- years[!(years %in% testyears[[j]])]
                        }
                        
                        crps <- matrix(nrow = 7, ncol = 7)
                        reslist <- list()
                        for(k in 1:7){
                          test <- testyears[[k]]
                          train <- c(1985:2012)[!c(1985:2012) %in% test]
                          
                          res <- evalmodels(station = "Axams",
                                            train = train,
                                            test = test,
                                            gamboost_cvr = TRUE)
                          
                          crps[k,] <- res$crps
                          reslist[[k]] <- res
                        }
                        
                        colnames(crps) <- names(res$crps)
                        return(reslist)
                      }
)

# extract CRPS
crps_cross <- matrix(nrow = nrep_cross, ncol = 7)
# loop over all repetitions
for(i in 1:length(res_cross)){
  #loop over all 7 folds (for 7 methods)
  crps_cross_int <- matrix(nrow = length(res_cross[[1]]), ncol = 7)
  for(j in 1:length(res_cross[[1]])){
    crps_cross_int[j,] <- res_cross[[i]][[j]]$crps
  }
  crps_cross[i,] <- colMeans(crps_cross_int, na.rm = TRUE)
}
colnames(crps_cross) <- names(res_cross[[1]][[1]]$crps) 

save(crps_cross, file = "crps_cross.rda")
save(res_cross, file = "res_cross.rda")
                
                
## boxplot of crps_cross for station Axams
boxplot(1 - crps_cross[,c(2,3,4)] / crps_cross[,6], ylim = c(-0.005, 0.065),
        names = c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS"),
        ylab = "CRPS skill score", col = "lightgray") 
abline(h = 0, col = pal["EMOS"], lwd = 2)
  
  



######################################################
### models learned on 24 years (1985-2008) and evaluated on 4 years (2009-2012)
### over all observation stations


data("StationsTyrol", package = "RainTyrol")
stations <- StationsTyrol$name
test <- 2009:2012
train <- 1985:2008


res_24to4_all <- applyfun(1:length(stations),
                          function(i){
                            
                            set.seed(7)
                            
                            res <- evalmodels(station = stations[i],
                                              train = train,
                                              test = test,
                                              gamboost_cvr = TRUE)
                            
                            return(res)
                          }
)

# extract crps
crps_24to4_all <- matrix(nrow = length(stations), ncol = 7)
# loop over all stations
for(i in 1:length(stations)){
  crps_24to4_all[i,] <- res_24to4_all[[i]]$crps
}

colnames(crps_24to4_all) <- names(res_24to4_all[[1]]$crps)
rownames(crps_24to4_all) <- stations

save(crps_24to4_all, file = "crps_24to4_all.rda")
save(res_24to4_all, file = "res_24to4_all.rda")



# boxplot of crps_24to4_all
s <- 1 - crps_24to4_all[, 2:4]/crps_24to4_all[,6] 
colnames(s) <- c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS")
matplot(t(s[,]), type = "l", lwd = 2, 
        col = gray(0.5, alpha = 0.2),
        lty = 1, axes = FALSE, 
        xlab = "", ylab = "CRPS skill score", xlim = c(0.5, 3.5))
lines(s[70,], col = "limegreen", type = "o", pch = 19, lwd = 2)  
# Station Axams is the 77th station which is the 70th in the list of complete stations
boxplot(s, add = TRUE, col = "transparent")
abline(h = 0, col = pal["EMOS"], lwd = 2)

   

## prepare data for map which shows where distforest performed better than gamlss or gamboostLSS based on the crps
# FIXME: data source
if(FALSE){  
  crps_map <- crps_24to4_all[,c("distforest", "gamlss", "gamboostLSS", "emos_log")]  
  
  # best method
  bst <- apply(crps_map, 1, which.min)
  
  # distance of forest to best other method
  dst <- crps_map[,1] - crps_map[cbind(1:nrow(crps_map), apply(crps_map[, -1], 1, which.min) + 1)]
  
  # breaks/groups
  brk <- c(-0.1, -0.05, -0.005, 0.005, 0.05, 0.1)
  #brk <- c(-0.1, -0.05, -0.01, 0.01, 0.05, 0.1)
  grp <- cut(dst, breaks = brk)
  
  # HCL colors (relatively flashy, essentially CARTO Tropic)
  clr <- colorspace::diverge_hcl(5, h = c(195, 325), c = 80, l = c(50, 90), power = 1.3)
  
  
  library("raster") # dem (digital elevation model)
  library("sp")     # gadm www.gadm.org/country
  
  
  data(StationsTyrol, package = "RainTyrol")
  data(MapTyrol, package = "RainTyrol")
  data(MapTyrol_border, package = "RainTyrol")
  # Create SpatialPointsDataFrame from station list
  sp <- SpatialPointsDataFrame(subset(StationsTyrol,
                                      select=c(lon,lat)),
                               data = subset(StationsTyrol,
                                             select = -c(lon,lat)),
                               proj4string = crs(dem))
  
  
  ## plot map of Tyrol with all 95 observations
  layout(cbind(1, 2), width = c(9, 1))
  par(mar = c(5,4,4,0.1))
  raster::image(dem, col = rev(gray.colors(100)),
                main="Stations in Tyrol", ylab = "Latitude", xlab = "Longitude", 
                xlim = c(9.8,13.2), 
                ylim = c(46.6, 47.87))
  plot(tirol, add = TRUE)
  points(sp, pch = c(21, 24, 25, 22)[bst], bg = clr[grp], col = "black", las = 1, cex = 1.5)
  legend(x = 9.8, y = 47.815, pch = c(21, 24, 25, 22), legend = c("Distributional forest", "Prespecified GAMLSS", "Boosted GAMLSS", "EMOS"), cex = 1, bty = "n")
  text(x = 10.3, y = 47.82, labels = "Models with lowest CRPS")
  mtext("CRPS\ndifference", side=4, las = TRUE, at = c(x = 13.5, y = 47.76), line = 0.3)
  par(mar = c(0.5,0.2,0.5,2.3))
  ## legend
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(-0.2, 0.2), xaxs = "i", yaxs = "i")
  rect(0, brk[-6], 0.5, brk[-1], col = rev(clr))
  axis(4, at = brk, las = 1, mgp=c(0,-0.5,-1))
}  
  