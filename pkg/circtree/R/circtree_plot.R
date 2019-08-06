plot.circtree <- function(x, terminal_panel = node_circular,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...){

  if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * 1
  if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
  partykit::plot.modelparty(x, terminal_panel = terminal_panel,
    tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
}


plot_circular <- function(X, coefs = NULL, stack = 10, cex = NULL, label = TRUE, 
  circlab = NULL, polygons = TRUE, rug = TRUE, 
  kernel_density = FALSE, type = c("response", "parameter"), 
  plot_type = c("mathematical", "geographics"), response_range = c(-pi, pi)){

  type <- match.arg(type)
  plot_type <- match.arg(plot_type)

  ## Convert to meteorological wind direction
  if(plot_type == "geographics"){
    X <- (-(X - pi/2)) %% (2 * pi)
    coefs[1] <- (-(coefs[1] - pi/2)) %% (2 * pi)
    circlab = c('E', 'N', 'W','S')
    if(is.null(cex)) cex <- 0.8
  }

  if((!is.null(circlab) | length(circlab) != 4) & plot_type != "geographics"){
    if(all(response_range == c(0, 2 * pi))){
      circlab <- c('$0$', '$\\pi/2$', '$\\pi$','$3/2\\pi$')
      if(is.null(cex)) cex <- 0.8
    } else if(all(response_range == c(0, 360))){
      circlab <- c('$0^\\degree$', '$90^\\degree$', '$180^\\degree$','$270^\\degree$')
      if(is.null(cex)) cex <- 0.6
    } else if(all(response_range == c(-pi, pi))){
      circlab <- c('$0$', '$\\pi/2$', '$\\pi$','$-\\pi/2$')
      if(is.null(cex)) cex <- 0.8
    } else if(all(response_range == c(-180, 360))){
      circlab <- c('$0^\\degree$', '$90^\\degree$', '$180^\\degree$','$-90^\\degree$')
      if(is.null(cex)) cex <- 0.6
    } else{
      circlab <- round(response_range, 2)
      if(is.null(cex)) cex <- 0.6
    }
  } else{
    if(is.null(cex)) cex <- 0.6
  }

  # Empty Plot
  par(mar = c(0, 0, 0, 0) + 0., cex = cex)
  #tmp <- par("pin")[1]/6
  #par(cex = tmp)
  plot(NA, xlim = c(-2,2), ylim = c(-2,2),
    type = "l", axes = FALSE, xlab = "",ylab = "", asp = 1)

  # Histogram
  Xhist <- hist(X, plot = FALSE, breaks = seq(0, 360, stack) * pi / 180)
  idx <- (Xhist$density != 0)
  Xscaled <- Xhist$density / max(Xhist$density)
 
  # Draw densities (lines or polygons)
  if (! polygons) {
    segments(cos(Xhist$mids[idx]), sin(Xhist$mids[idx]),
             (Xscaled[idx] + 1) * cos(Xhist$mids[idx]),
             (Xscaled[idx] + 1) * sin(Xhist$mids[idx]))
    #points(1 * cos(Xhist$mids), 1 * sin(Xhist$mids), pch = 20, col = c(NA, gray(0.2, alpha = 0.7))[idx + 1]) 
  } else {
    segwidth <- diff(Xhist$mids[1:2])     # width of the segment
    segn     <- ceiling(segwidth / 0.005) # level of detail for the curvature of the segment
    for (i in which(idx)) {
      x <- seq(Xhist$mids[i] - segwidth/2, Xhist$mids[i] + segwidth/2, length = segn)
      polygon(c(cos(x), (Xscaled[i] + 1) * cos(rev(x))),
              c(sin(x), (Xscaled[i] + 1) * sin(rev(x))), border = NA, col = "gray80")
    }
  } 

  # Plot rugs
  if(rug){
    segments(cos(X), sin(X),
             (0.9) * cos(X),
             (0.9) * sin(X), gray(0.2, alpha = 0.7))
  }

  # Plot circle
  circ_ln <- seq(0, 360, 0.5) * pi / 180
  lines(cos(circ_ln), sin(circ_ln))
  
  # Plot ticks and labels
  circ_sgm <- seq(0, 360, 45) * pi / 180
  segments(.8 * cos(circ_sgm), .8 * sin(circ_sgm), 1 * cos(circ_sgm), 1 * sin(circ_sgm))

  if(label & requireNamespace("latex2exp", quietly = TRUE)){
    text(.8, 0, latex2exp::TeX(circlab[1]), adj = c(1, 0.5))
    text(0, .8, latex2exp::TeX(circlab[2]), adj = c(0.5, 1))
    text(-.8, 0, latex2exp::TeX(circlab[3]), adj = c(0, 0.5))
    text(0, -.8, latex2exp::TeX(circlab[4]), adj = c(0.5, 0))
  }


  # Plot kernel density
  if(kernel_density){  
    d <- circular::density.circular(circular::circular(X), bw = 1)
    lines(d, col = 2, lty = 2, lwd = 1.5)
  }

  # Plot fitted density
  if(!is.null(coefs)){
    arrows(0, 0, 0.5 * cos(coefs[1]), 0.5 * sin(coefs[1]), col = 2, length = 0, lwd = 2)
    #arrows(0, 0, 0.5 * cos(coefs[1]), 0.5 * sin(coefs[1]), col = 2, length = 0.1, code = 2, lwd = 2)
    circular::plot.function.circular(function(x) 
      dvonmises(x, circular::circular(coefs[1]), coefs[2]), add = TRUE, col = 2, lty = 1, lwd = 1.5)
  }
}


node_circular <- function(obj, which = NULL, id = TRUE, pop = TRUE,
  xlab = FALSE, ylab = FALSE, mainlab = NULL, type = c("response", "parameter"), 
  plot_type = c("mathematical", "geographics"), ...){

  type <- match.arg(type)
  plot_type <- match.arg(plot_type)

  ## obtain dependent variable on range of (0, 2*pi)
  response_range <- attr(obj$fitted[["(response)"]], "response_range")
  y <- angle_retrans(obj$fitted[["(response)"]], 0, 2*pi)
  fitted <- obj$fitted[["(fitted)"]]

  ## y- and x-lab
  if(isTRUE(ylab)) ylab <- names(y)
  if(identical(ylab, FALSE)) ylab <- ""
  if(identical(xlab, FALSE)) xlab <- ""

  ## set up appropriate panel functions
  rval <- function(node) {
    
    ## node index
    nid <- partykit::id_node(node)
    ix <- fitted %in% partykit::nodeids(obj, from = nid, terminal = TRUE)

    ## dependent variable
    y <- y[ix]
    if(length(y) > 1000) y <- sample(y, 1000)

    ## coefficients
    coefs <- c(signif(coef(node$info$object)[1],2), signif(coef(node$info$object)[2],2))

    ## set up viewport tree
    top.vp <- grid::viewport(layout = grid::grid.layout(ncol = 1, nrow = 3,
             widths = grid::unit(c(1, 1, 1), c("null", "null", "null")), heights = grid::unit(c(0.3, 0.1, 0.6), c("npc", "npc", "npc"))),
             width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), name = paste("node", nid, sep = ""))
    topmargin <- grid::viewport(layout.pos.col = 1, layout.pos.row = 1, name = "topmargin")
    text <- grid::viewport(layout.pos.col = 1, layout.pos.row = 2, name = "text")
    plot <- grid::viewport(layout.pos.col = 1, layout.pos.row = 3, name = "plot", )
    
    splot <- grid::vpTree(top.vp, grid::vpList(topmargin, text, plot))
    grid::pushViewport(splot)

    ## main title
    grid::seekViewport("text")

    if (is.null(mainlab)) { 
      mainlab <- if(id) {
      function(id, nobs, mu, kappa, type){
        ## Convert to response range if type = 'response'
        if(type == "response"){
          mu <- signif(angle_retrans(mu, response_range[1], response_range[2]), 2)
        }
        sprintf("Node %s (n = %s) \n mu = %s, kappa = %s", id, nobs, mu, kappa)
      }
      #function(id, nobs, mu, kappa) sprintf("Node %s (n = %s)", id, nobs)
      } else {
      function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, node$info$object$ny, coefs[1], coefs[2], type)
    }
    grid::grid.text(mainlab, y = grid::unit(0.9, "npc"))

    ## plot rectangle and actual graphic, optional x and ylab 
    grid::seekViewport("plot")

    grid::grid.rect(gp = grid::gpar(fill = "transparent", col = 1), width = grid::unit(0.9, "npc"))
    gridGraphics::grid.echo(function() plot_circular(y, coefs, plot_type = plot_type, 
      response_range = response_range, ...), newpage = FALSE)

    if(ylab != "") grid::grid.text(ylab, y = grid::unit(0.5, "npc"), x = grid::unit(-2.5, "lines"), rot = 90)
    if(xlab != "") grid::grid.text(xlab, x = grid::unit(0.5, "npc"), y = grid::unit(-2, "lines"))         

    if(pop) grid::popViewport() else grid::upViewport()
    if(pop) grid::popViewport() else grid::upViewport()
  }
  
  return(rval)
}
class(node_circular) <- "grapcon_generator"

