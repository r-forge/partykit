
library("partykit")
set.seed(290875)

### regression
airq <- subset(airquality, complete.cases(airquality))
airct <- ctree(Ozone ~ ., data = airq)

mean((airq$Ozone - predict(airct))^2)
#logLik(airct)
#logLik(airct, airq, perm = "Temp")

varimp(airct)

aircf <- cforest(Ozone ~ ., data = airq, ntree = 100)

varimp(aircf)

varimp(aircf, conditional = TRUE)

ict <- cforest(Species ~ ., data = iris, ntree = 100)
varimp(ict)
varimp(ict, risk = "misclass")

set.seed(29)
varimp(ict, risk = "misclass", conditional = TRUE)

