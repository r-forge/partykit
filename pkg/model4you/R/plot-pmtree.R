
#' Survival plot for a given survreg model
#'
#' Can be used on its own but is also useable as plotfun in 
#' \code{\link{node_pmterminal}}.
#'
#' @param mod A model of class survreg.
#' @param data optional data frame. If NULL the data stored in mod is used.
#' @param theme A ggplot2 theme.
#'
#' @return RETURN DESCRIPTION
#' @examples
#' if(require("survival")) {
#'   survplot(survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian))
#' }
#' 
#' @importFrom ggplot2 ggplot geom_line theme_classic aes_string
#' @importFrom stats model.frame predict formula terms get_all_vars
#' @importFrom survival Surv survreg
#' @importFrom Formula as.Formula 
#' @export
survplot <- function(mod, data = NULL, theme = theme_classic()) {
  
  ## get formula and data
  modcall <- getCall(mod)
  modformula <- as.Formula(eval(modcall$formula))
  xformula <- formula(modformula, lhs = 0, rhs = 1)
  if(is.null(data)) data <- eval(modcall$data)
  xdat <- unique(get_all_vars(xformula, data = data))
  
  ## get survivor functions for each treatment group
  p <- seq(.01, .99, by=.01)
  pr_raw <- predict(mod, newdata = xdat,
                    type = "quantile", p = p)
  pr <- do.call("rbind", 
                lapply(1:NROW(xdat), 
                       function(i) data.frame(xdat[i, , drop = FALSE],
                                              pr = pr_raw[i, ], prob = rev(p),
                                              row.names = NULL)))
  
  ## plot
  xnam <- attr(terms(xformula), "term.labels")
  ggplot(data = pr, aes_string(x = "pr", y = "prob", group = xnam, 
                               color = xnam)) + 
    geom_line() + theme
}



#' Panel-Generator for Visualization of pmtrees
#'
#' The plot method for party and constparty objects are rather flexible and can 
#' be extended by panel functions. The pre-defined panel-generating function of 
#' class grapcon_generator for pmtrees is documented here. 
#'
#' @param obj an object of class party.
#' @param digits integer, used for formating numbers.
#' @param confint Should a confidence interval be computed.
#' @param plotfun Plotting function to be used. Needs to be of format 
#' \code{function(mod, data)} where \code{mod} is the model object. 
#' See examples for more details.
#' @param nid function to retrieve info on what is plottet as node ids.
#'
#' @examples
#' if(require("survival")) {
#' ## compute survreg model
#' mod_surv <- survreg(Surv(futime, fustat) ~ factor(rx), ovarian, dist='weibull')
#' survplot(mod_surv)
#' 
#' ## partition model and plot
#' tr_surv <- pmtree(mod_surv)
#' plot(tr_surv, terminal_panel = node_pmterminal(tr_surv, plotfun = survplot, 
#'                                                confint = TRUE))
#' }
#' 
#' @export
#' 
#' @importFrom grid viewport grid.layout unit popViewport 
#' pushViewport grid.rect grid.draw gpar grid.text upViewport
#' @importFrom gridExtra tableGrob ttheme_minimal 
#' @importFrom partykit id_node
node_pmterminal <- function(obj, digits = 2, confint = TRUE, plotfun, 
                            nid = function(node) paste0(nam[id_node(node)], ", n = ", node$info$nobs))
{
  
  nam <- names(obj)
  mod <- obj$info$object
  wterminals <- predict(obj, type = "node")
  dat <- obj$data
  
  ### panel function for the inner nodes
  rval <- function(node, .nid = nid) {  
    
    nid <- .nid(node)
    
    ## model
    nmod <- update(mod, data = dat, subset = (wterminals == id_node(node)))
    coefs <- as.matrix(node$info$coefficients)
    if(confint) {
      ci <- confint(nmod)
      coefs <- cbind(coefs, ci)
    }
    cf <- format(round(coefs, digits), nsmall = digits)
    colnams <- colnames(cf)
    colnams[1] <- "theta"
    cftab <- tableGrob(cf, cols = colnams,
                       theme = ttheme_minimal(colhead = list(fg_params = list(parse=TRUE))))
    tabwid <- sum(cftab$widths)
    
    node_vp <- viewport(
      layout = grid.layout(2, 1),
      width = max(tabwid, unit(0.95, "npc")),
      height = unit(0.95, "npc")
    )
    pushViewport(node_vp)
    
    grid.rect()
    
    ## table
    tablevp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
    pushViewport(tablevp)
    grid.draw(cftab)
    popViewport()
    
    ## plot
    plotvp <- viewport(layout.pos.row = 2, layout.pos.col = 1,
                       width = unit(0.95, "npc"),
                       height = unit(0.95, "npc"))
    pushViewport(plotvp)
    pl <- plotfun(nmod, data = dat)
    print(pl, vp = plotvp)
    popViewport()
    
    
    # nidtxt <- paste0(nam[nid], ", n = ", node$info$nobs)
    nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
                         width = max(unit(1, "lines"), unit(1.3, "strwidth", nid)),
                         height = max(unit(1, "lines"), unit(1.3, "strheight", nid)))
    pushViewport(nodeIDvp)
    grid.rect(gp = gpar(fill = "white"))
    grid.text(nid)
    popViewport()
    
    upViewport()
  }
  
  return(rval)
}

class(node_pmterminal) <- "grapcon_generator"