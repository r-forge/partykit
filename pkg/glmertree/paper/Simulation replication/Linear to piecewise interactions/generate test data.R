N <- list() # different levels for N
N[[1]] <- 200
N[[2]] <- 500
N[[3]]<- 1000

bi <- list() # different numbers and values for the study intercepts
bi$numbclus <- c(25)
bi$sigmas <- c(2.5, 7.5)

rho <- list() # levels of intercorrelations between the covariates
rho[[1]] <- matrix(0, nrow=20, ncol=20)
rho[[2]] <- matrix(.3, nrow=20, ncol=20)
diag(rho[[1]]) <- 1
diag(rho[[2]]) <- 1

corUbi <- list()
corUbi[[1]] <- "uncorrelated"
corUbi[[2]] <- "bi and U correlated"

pXj <- list() # levels of number of covariates
pXj[[1]] <- 5; pXj[[2]] <- 10 

types <- list("linear","both","piecewise") # different levels for the difference in treatment outcome


testdata <- list()
testdescriptions <- list()

set.seed(38916032)

for (c in 1:50) {  
  counter <- 0
  for (d in 1:length(bi$sigmas)) {# d is counter for sigma_b value
    for (e in 1:length(bi$numbclus)) {# e is counter for number of clusters
      numbclus <- bi$numbclus[[e]] # numbclus is number of different random intercept values
      for (f in 1:length(corUbi)) { # f is counter for correlation between U and b
        for (g in 1:length(types)) { # g is counter for type
          type <- types[[g]]
          for (h in 1:length(pXj)) { # h is counter for number of covariates
            p <- pXj[[h]]
            for (i in 1:length(N)) { # i is counter for sample size N
              nobs <- N[[i]]
              for (j in 1:length(rho)) { # j is counter for intercorrelations between covariates
                L <- chol(rho[[j]]) # Cholesky decomposition of correlation matrices
                nvars <- dim(rho[[j]])[1]
                
                # Random variables that follow an M correlation matrix
                r <- t(L) %*% matrix(rnorm(nvars*nobs, sd=10), nrow=nvars, ncol=nobs)
                r <- t(r)
                
                # create centered covariate datasets
                r <- as.data.frame(r)
                names(r) <- paste("X", 1:15, sep="")
                
                # add population means to covariates
                r$X1 <- r$X1+10 # population mean of covariate X1
                r$X2 <- r$X2+30 # population mean of covariate X2
                mu3 <- round(runif(1, -70, 70))
                r$X3 <- r$X3+mu3# population mean of covariate X3
                r$X4 <- r$X4-40 # population mean of covariate X4
                r$X5 <- r$X5+70 # population mean of covariate X5
                mu6t15 <- ceiling(runif(15, -71, 70)) # population mean of covariate X6 through X20
                for (mu in 6:15) {
                  r[,paste("X",mu,sep="")] <- r[,paste("X",mu,sep="")] + mu6t15[mu-5]
                }
                
                # add treatment variable to datasets
                r$T <- rbinom(N[[i]], 1, .5)
                
                if(type == "linear" | type == "both") {
                  # Main effect of X2:
                  r$Y <- -45 + 1.5*r$X2
                  # Centered interactions:
                  r$Y <- r$Y + (r$X2-30)*(r$X1-10)*-.25
                  r$Y <- r$Y + (r$X2-30)*(r$X5-70)*.25
                  # Centered treatment interactions:
                  r$Y <- r$Y + r$T*(r$X2-30)*(r$X1-10)*-.25
                  r$Y <- r$Y + r$T*(r$X2-30)*(r$X5-70)*.25
                } else {
                  node3 = r$X2<=30 & r$X1<=17
                  node4 = r$X2<=30 & r$X1>17
                  node6 = r$X2>30 & r$X5<=63
                  node7 = r$X2>30 & r$X5>63
                  r$Y[node3 & r$T == 0] <- -20.17224
                  r$Y[node3 & r$T == 1] <- -28.39720
                  r$Y[node4 & r$T == 0] <- 13.75540
                  r$Y[node4 & r$T == 1] <- 39.59950
                  r$Y[node6 & r$T == 0] <- -13.75933
                  r$Y[node6 & r$T == 1] <- -39.53344
                  r$Y[node7 & r$T == 0] <- 20.14495
                  r$Y[node7 & r$T == 1] <- 28.37829
                }
                if(type == "both") {
                  r$Y <- .5*r$Y
                  node3 = r$X2<=30 & r$X1<=17
                  node4 = r$X2<=30 & r$X1>17
                  node6 = r$X2>30 & r$X5<=63
                  node7 = r$X2>30 & r$X5>63
                  r$Y[node3 & r$T == 0] <- r$Y[node3 & r$T == 0] + .5*-20.17224
                  r$Y[node3 & r$T == 1] <- r$Y[node3 & r$T == 1] + .5*-28.39720
                  r$Y[node4 & r$T == 0] <- r$Y[node4 & r$T == 0] + .5*13.75540
                  r$Y[node4 & r$T == 1] <- r$Y[node4 & r$T == 1] + .5*39.59950
                  r$Y[node6 & r$T == 0] <- r$Y[node6 & r$T == 0] + .5*-13.75933
                  r$Y[node6 & r$T == 1] <- r$Y[node6 & r$T == 1] + .5*-39.53344
                  r$Y[node7 & r$T == 0] <- r$Y[node7 & r$T == 0] + .5*20.14495
                  r$Y[node7 & r$T == 1] <- r$Y[node7 & r$T == 1] + .5*28.37829
                }
                
                
   
                # add study intercept variable to datasets, which is randomly correlated with X1 thru X5
                sigma_bi <- bi$sigmas[d]
                #if(sigma_bi==0) {int_vals <- rep(rnorm(numbclus,0,300),each=nobs/numbclus)}
                if(sigma_bi!=0) {int_vals <- rep(rnorm(numbclus,0,bi$sigmas[d]),each=nobs/numbclus)} # generate values to be use as random intercepts
                int_vals <- int_vals[order(int_vals)] # order random intercept values
                
                if (corUbi[[f]]=="bi and U correlated") {
                  correlated_covariate <- sample(1:5, 1)
                  r$bi <- r[,paste("X", correlated_covariate, sep="")] + rnorm(nobs, sd=20)
                }
                if (corUbi[[f]]=="uncorrelated") {
                  correlated_covariate <- NA
                  r$bi <- rnorm(nobs)
                }
                r$bi[order(r$bi)] <- int_vals
                if(sigma_bi==0){r$Y <- r$Y}
                if(sigma_bi!=0){r$Y <- r$Y + r$bi}
                r$cluster <- factor(r$bi)
                levels(r$cluster) <- 1:length(unique(r$bi))
                
                # generate and add error to Y
                r$errorY <- rnorm(N[[i]], sd=5)
                r$Y <- r$Y + r$errorY
                
                
              # write dataset to list and add description
              counter <- counter+1
              testdata[[counter]] <- list() 
              r$T <- factor(r$T)
              testdata[[counter]][[3]] <- r[,c(paste("X",1:p, sep=""), "T", "bi", "cluster", "errorY", "Y")]
              r$T <- 0 # dataset[[1]] has T=1
              r$T <- factor(r$T)
              testdata[[counter]][[1]] <- r[,c(paste("X",1:p, sep=""), "T", "cluster", "bi")]
              r$T <- 1 # dataset[[2]] has T=2
              r$T <- factor(r$T)
              testdata[[counter]][[2]] <- r[,c(paste("X",1:p, sep=""), "T", "cluster", "bi")]
              testdescriptions[[counter]] <- list(
                paste("N =", N[[i]]),
                paste("rho =", rho[[j]][1,2]), 
                paste("number of covariates =", pXj[[h]]),
                paste("type =", type),
                paste("correlation between U and bi =", corUbi[[f]]),
                paste("no. of random intercept values =", numbclus),
                paste("sigma bi =", bi$sigmas[d]),
                paste("random intercept is correlated with X", correlated_covariate, sep=""))
             }
            }
          }
        }
      }
    }
  }
  save(testdata, file=paste("testdatasets", c, sep=""))
  save(testdescriptions, file=paste("testdescriptions", c, sep=""))
}