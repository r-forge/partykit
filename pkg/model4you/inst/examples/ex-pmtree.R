if(require("TH.data") & require("survival")) {
  ## base model
  bmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2, model = TRUE)
  survreg_plot(bmod)
  
  ## partitioned model
  tr <- pmtree(bmod)
  plot(tr, terminal_panel = node_pmterminal(tr, plotfun = survreg_plot, 
                                            confint = TRUE))
  summary(tr)
  summary(tr, node = 1:2)
  
  logLik(bmod)
  logLik(tr)
}
