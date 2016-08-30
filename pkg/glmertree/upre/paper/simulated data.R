###################
## Simulate data ##
###################

set.seed(42)
simdata <- data.frame(
  x1 = runif(500, -10, 30),
  x2 = runif(500, 35, 50),
  x3 = runif(500, 10, 12),
  x4 = rnorm(500, 5, 5),
  x5 = rnorm(500, 0, 10),
  x6 = rnorm(500, -10, 15),
  x7 = rnorm(500, sd = 7.5),
  x8 = rnorm(500, sd = 7.5),
  x9 = runif(500, 8, 14),
  x10 = rnorm(500, -5, 7)
)
rules <- list(
  rule1 = "x1 > 5 & x2 < 40",
  rule2 = "x3 > 11.25 & x4 < 6",
  rule3 = "x5 > 2 & x6 < -12.5"
)
error <- rnorm(500, 0, 5)
for(i in 1:length(rules)) {
  simdata[[paste("rule",i,sep="")]] <- as.numeric(with(simdata, eval(parse(text = rules[[i]]))))
}

simdata$y <- with(simdata, x8*.15 + x9*.8 + rule1*15 + rule2*15 + rule3*15) + error
cor(simdata$y, with(simdata, x8*.15 + x9*.8 + rule1*3 + rule2*3 + rule3*3))
cor(simdata$y, error)
apply(simdata, 2, sd)
cor(simdata[,c(paste("rule",1:3,sep=""),paste("x",7:10,sep=""),"y")])
simdata2 <- simdata[,c(paste("x", 1:10, sep=""), "y")]
simdata2 <- round(simdata2, digits = 1)