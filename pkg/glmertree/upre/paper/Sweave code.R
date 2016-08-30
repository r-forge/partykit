setwd("C:/Users/tobii/Desktop/swReg/upre/paper")
Sweave("prediction_rule_ensembles.Rnw")
tools::texi2dvi("prediction_rule_ensembles.tex", pdf = TRUE, clean = TRUE)
