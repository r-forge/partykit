if(require("TH.data") & require("survival")) {
  ## base model
  bmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2, model = TRUE)
  survplot(bmod)
  
  ## partitioned model
  tr <- pmtree(bmod)
  plot(tr, terminal_panel = node_pmterminal(tr, plotfun = survplot, 
                                            confint = TRUE))
  summary(tr)
  summary(tr, node = 1:2)
  
  logLik(bmod)
  logLik(tr)
}
