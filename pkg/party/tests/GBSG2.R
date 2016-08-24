
library("partyNG")
library("survival")

data("GBSG2", package = "TH.data")

ctree(Surv(time, cens) ~ horTh + age + menostat + 
    tsize + tgrade + pnodes + progrec + estrec, data = GBSG2)

ctree(Surv(time, cens) ~ .,  data = GBSG2)

partykit::ctree(Surv(time, cens) ~ horTh + age + menostat +
    tsize + tgrade + pnodes + progrec + estrec, data = GBSG2)
