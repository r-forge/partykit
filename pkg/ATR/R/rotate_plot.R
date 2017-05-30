
rotate <- function(m, to = "left", ...) {
  stopifnot(to == "left")
  if(inherits(m, c('constparty', 'party'))) {
    class(m) <- c('left_rotated_tree', class(m))
    m
  }else{
    warning("Could not typecast to ", sQuote("left_rotated_tree"), 
            " ! Returning unmodified object")
    m
  }
}

.plot_node.left_rotated_tree <- function(node, obj, xlim, ylim, nx, ny, terminal_panel, inner_panel, edge_panel,
                                         tnex=2, drop_terminal=TRUE, ...)
{

  ### unpacking '...'
  kwargs <- list(...)
  if(!is.null(kwargs$debug)) {debug <- kwargs$debug} else {debug <- FALSE}
  if(!is.null(kwargs$cex)) {cex <- kwargs$cex} else {cex <- 1} # removes all x-axis except of the last one
  if(!is.null(kwargs$remove.xaxis)) {remove.xaxis <- kwargs$remove.xaxis} else {remove.xaxis <- FALSE}
  if(!is.null(kwargs$tree.offset)) {tree.offset <- kwargs$tree.offset} 
  else {tree.offset <- 0} #shifts the entire tree x lines up
  if(!is.null(kwargs$remove.nobs)) {remove.nobs <- kwargs$remove.nobs} else {remove.nobs <- FALSE}
  if(!is.null(kwargs$nobs.loc)) {nobs.loc <- kwargs$nobs.loc} else {nobs.loc <- 'right'} #nobs location
  #options
  
  
  ### the workhorse for plotting trees
  
  ### set up viewport for terminal node
  if (is.terminal(node)) {
    y <- ylim[2]
    x <- xlim[2] - 0.5*tnex # starting from right plot is centered at half of its width
    
    if(remove.xaxis && ylim[1]<1) { #This if clause is buggy
        tn_vp <- viewport(y = unit(y, "native"),
                          x = unit(x, "native"),
                          height = unit(1, "native")+unit(2*cex,"lines"), 
                          width =  unit(tnex, "native"),
                          just = c("center", "top"),# just = c("left", "center")
                          name = paste("Node", id_node(node), sep = ""),
                          gp=gpar(cex=cex))
        pushViewport(tn_vp)
        if (debug) grid.rect(gp = gpar(fill=FALSE, lty = "dotted", col = 4))
        terminal_panel(node) 
      } else {
        tn_vp <- viewport(y = unit(y, "native"),
                          x = unit(x, "native"),
                          height = unit(1, "native"), 
                          width =  unit(tnex, "native"),
                          just = c("center", "top"),# just = c("left", "center")
                          name = paste("Node", id_node(node), sep = ""),
                          gp=gpar(cex=cex))
        pushViewport(tn_vp)
        if (debug) grid.rect(gp = gpar(lty = "dotted", col = 4))
        terminal_panel(node) 
      }

    upViewport()
    return(NULL)
  }    
  
  ## convenience function for computing relative position of splitting node
  pos_frac <- function(node) {
    if(is.terminal(node)) 0.5 else {
      width_kids <- sapply(kids_node(node), width)
      nk <- length(width_kids)
      rval <- if(nk %% 2 == 0) sum(width_kids[1:(nk/2)]) else
        mean(cumsum(width_kids)[nk/2 + c(-0.5, 0.5)])
      rval/sum(width_kids)
    }
  }
  
  ## extract information
  split <- split_node(node)
  kids <- kids_node(node)
  width_kids <- sapply(kids, width)
  nk <- length(width_kids)
  
  ### position of inner node
  y0 <- ylim[1] + pos_frac(node) * diff(ylim)
  x0 <- min(xlim)
  
  ### relative positions of kids
  yfrac <- sapply(kids, pos_frac)
  y1lim <- ylim[1] + cumsum(c(0, width_kids))/sum(width_kids) * diff(ylim)
  y1 <- y1lim[1:nk] + yfrac * diff(y1lim)
  if (!drop_terminal) {
    x1 <- rep(x0 + 1, nk)
  } else {
    x1 <- ifelse(sapply(kids, is.terminal), xlim[2] - tnex, x0 + 1)
  }
  
  ### draw edges
  grid.lines(x = unit(c(x0, x0+.5), "native"), 
             y = unit(c(y0, y0), "native") + unit(tree.offset, "native")) #horizontal
  
  terminal.kid = sapply(kids, is.terminal) #boolean list containing whether child is terminal or not.
  for(i in 1:nk) {
    grid.lines(x = unit(c(x0+.5, x0+.5), "native"), 
               y = unit(c(y0, y1[i]), "native") + unit(tree.offset, "native")) #vertical
    
    if(terminal.kid[i]) {
      
      #if (debug) grid.points(x=unit(x1[i]-0.5, "native"), 
      #                       y=unit(y1[i]+0.5, "native") - unit(1, "lines"),
      #                       pch=4)
      grid.lines(x = unit(c(x0+.5, xlim[2]-tnex), "native")-unit(c(0,.4),"lines"), 
                 y = unit(c(y1[i], y1[i]), "native") + unit(tree.offset, "native")) #horizontal
      # grid.lines(x = unit(c(x0+.5, x1[i]-0.5), "native"), 
      #            y = unit(c(y1[i]+.5, y1[i]+.5), "native")-unit(1, "lines")) #horizontal
    } else {
      grid.lines(x = unit(c(x0+.5, x1[i]), "native"), 
                 y = unit(c(y1[i], y1[i]), "native") + unit(tree.offset, "native")) #horizontal
    }
  } 
  
  ### position of labels
  xpos <- x0 + 0.5
  ypos <- y0 - (y0 - y1)/2
  
  ### setup labels
  for(i in 1:nk) {
    if(terminal.kid[i]) {
      w = unit(xlim[2]-tnex-xpos, "native") +unit(0.6*cex, "lines")
    } else {
      w = unit(1, "native") - unit(.4*cex, "lines")
    }
    sp_vp <- viewport(x = unit(xpos, "native")-unit(1*cex, "lines"),
                      y = unit(ypos[i], "native")+unit(.25*(-1)^(i+1)*cex,"lines") + unit(tree.offset, "native"),
                      width = w,
                      height = unit(1*cex, "lines"), just= "left",
                      name =  paste("edge", id_node(node), "-", i, sep = ""),
                      gp=gpar(cex=cex))
    pushViewport(sp_vp)
    if(debug) grid.rect(gp = gpar(fill=FALSE, lty = "dotted", col = 2))
    edge_panel(node, i)
    upViewport()
  }
  
  ### node ids for all nodes
  fill <- "white"
  fill <- rep(fill, length.out = 2L)
  for(i in 1:length(kids)) {
    nodeID <- viewport(x = unit(x0 + .5, "native"), # - unit(.5, "lines")
                       y = unit(y1[i], "native") + unit(tree.offset, "native"),
                       width = max(unit(1*cex, "lines"), unit(1.3*cex, "strwidth", as.character(id_node(kids[[i]])))),
                       height = max(unit(1*cex, "lines"), unit(1.3*cex, "strheight", as.character(id_node(kids[[i]])))),
                       just="right", gp=gpar(cex=cex)
    )
    pushViewport(nodeID)
    grid.rect(gp = gpar(fill = fill[2]), just = "center")
    grid.text(as.character(id_node(kids[[i]])), just = "center")
    popViewport()
  }
  
  ### number of observations
  if(!remove.nobs) {
    for(i in 1:nk) {
      #if(!terminal.kid[i]) {
        ## extract data
        nid <- id_node(kids[[i]])
        dat <- data_party(obj, nid)
        yn <- dat[["(response)"]]
        wn <- dat[["(weights)"]]
        n.text <- sprintf("n = %s", sum(wn))
        if(is.null(wn)) wn <- rep(1, NROW(yn))
        if(nobs.loc == 'top') {
          #possibility to display number of observation above/below node number
          nobs_vp <- viewport(x = unit(x0 + .5, "native") - max(unit(.5*cex, "lines"), 
                                                                unit(.65*cex, "strwidth", as.character(nid))),
                              y = unit(y1[i], "native") + unit(.5*(-1)^i*cex, "lines")
                                  + (-1)^i*max(unit(.5*cex, "lines"), 
                                        unit(.65*cex, "strheight", as.character(nid)))
                                  + unit(tree.offset, "native"),
                              width = unit(1*cex, "strwidth", n.text),
                              height = unit(1*cex, "lines"), just ="center",
                              name =  paste("nobs", id_node(node), "-", i, sep = ""),
                              gp=gpar(cex=cex))
        } else {
          #if nobs.loc == 'right'
          nobs_vp <- viewport(x = unit(x0 + .5, "native") + unit(.2*cex, "lines"),
                              y = unit(y1[i], "native") + unit(.5*(-1)^i*cex, "lines") 
                                  + unit(tree.offset, "native"),
                              width = unit(1*cex, "strwidth", n.text),
                              height = unit(1*cex, "lines"), just ="left",
                              name =  paste("nobs", id_node(node), "-", i, sep = ""),
                              gp=gpar(cex=cex))
        }
        pushViewport(nobs_vp)
      
        if(debug) grid.rect(gp = gpar(fill=FALSE, lty = "dotted", col = 2))
        grid.text(n.text)
        upViewport()
      #}
    }
  }
  
  ### create viewport for inner node
  in_vp <- viewport(x = unit(x0, "native"),
                    y = unit(y0, "native") + unit(tree.offset, "native"),
                    height = unit(.5, "native"),
                    width = unit(.3, "native"), 
                    name = paste("Node", id_node(node), sep = ""),
                    gp=gpar(cex=cex))
  pushViewport(in_vp)
  if(debug) grid.rect(gp = gpar(fill=FALSE, lty = "dotted"))
  inner_panel(node, ...)
  upViewport()
  
  ## call workhorse for kids
  for(i in 1:nk) 
    .plot_node.left_rotated_tree(kids[[i]], obj,
                               c(x1[i], xlim[2]), c(y1lim[i], y1lim[i+1]), nx, ny, 
                               # Note: this might be wrong! I'd expect something else than nx
                               terminal_panel, inner_panel, edge_panel,
                               tnex = tnex, drop_terminal = drop_terminal, ...)
}

plot.left_rotated_tree <- function(x, main = NULL,
                                 terminal_panel = node_ecdf, tp_args = list(),
                                 inner_panel = node_inner.left_rotated_tree, ip_args = list(),
                                 edge_panel = edge_simple.left_rotated_tree, ep_args = list(),
                                 type = "extended", drop_terminal = NULL, tnex = NULL, 
                                 newpage = TRUE, pop = TRUE, gp = gpar(), ...)
{

  # unpacking '...'
  kwargs <- list(...)
  if(!is.null(kwargs$debug)) {debug <- kwargs$debug} else {debug <- "FALSE"}
  
  obj <- x
  ### compute default settings
  type <- match.arg(type)
  stopifnot(type == "extended")
  if (is.null(tnex)) tnex <- 2
  if (is.null(drop_terminal)) drop_terminal <- TRUE
  
  ### extract tree
  node <- node_party(x)
  ### total number of terminal nodes
  ny <- width(node)
  ### maximal depth of the tree
  nx <- depth(node, root = TRUE)
  
  ## setup newpage
  if (newpage) grid.newpage()
  
  ## setup root viewport
  root_vp <- viewport(layout = grid.layout(3, 3, 
                                           heights = unit(c(ifelse(is.null(main), 0, 2), 1, 2), 
                                                          c("lines", "null", "lines")),
                                           widths = unit(c(1, 1, 1), 
                                                         c("lines", "null", "lines"))), 
                      name = "root",
                      gp = gp)       
  pushViewport(root_vp)
  
  ## viewport for main title (if any)
  if (!is.null(main)) {
    main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1, 
                        name = "main")
    pushViewport(main_vp)
    if(debug) grid.rect(gp = gpar(fill=FALSE, lty = "dotted", col = 3))
    grid.text(y=unit(1, "lines"), main, just = "center")
    upViewport()
  }
  
  ## setup viewport for tree
  tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
                      yscale = c(0, ny), xscale = c(0, nx + (tnex - 1)), 
                      name = "tree")
  pushViewport(tree_vp)
  if(debug) grid.rect(gp = gpar(fill=FALSE, lty = "dotted", col = 3))
  
  ### setup panel functions (if necessary)
  if(inherits(terminal_panel, "grapcon_generator"))
    terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
  if(inherits(inner_panel, "grapcon_generator"))
    inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
  if(inherits(edge_panel, "grapcon_generator"))
    edge_panel <- do.call("edge_panel", c(list(x), as.list(ep_args)))
  
  
  if((nx <= 1 & ny <= 1)) {
    pushViewport(plotViewport(margins = rep(1.5, 4), name = paste("Node", id_node(node), sep = "")))
    terminal_panel(node)
  } else {
    
    ## call the workhorse
    .plot_node.left_rotated_tree(node, obj,
                               ylim = c(0, ny), xlim = c(0.25, nx + (tnex - 1)),
                               nx = nx, ny = ny, 
                               terminal_panel = terminal_panel,
                               inner_panel = inner_panel,
                               edge_panel = edge_panel,
                               tnex = tnex,
                               drop_terminal = drop_terminal, ...)
  }
  upViewport()
  if (pop) popViewport() else upViewport()
}

### Plot function for inner node ###
# To change it, write a new function (change class to "grapcon_generator") and execute
# R> plot.left_rotated_tree(..., inner_panel = YOURFUNCTION)
node_inner.left_rotated_tree <- function(obj, id = TRUE, pval = TRUE, abbreviate = FALSE, fill = "white", gp = gpar())
{

  meta <- obj$data
  nam <- names(obj)
  
  extract_label <- function(node) {
    if(is.terminal(node)) return(rep.int("", 2L))
    
    varlab <- character_split(split_node(node), meta)$name
    if(abbreviate > 0L) varlab <- abbreviate(varlab, as.integer(abbreviate))
    
    ## FIXME: make more flexible rather than special-casing p-value
    if(pval) {
      pval <- suppressWarnings(try(!is.null(info_node(node)$p.value), silent = TRUE))
      pval <- if(inherits(pval, "try-error")) FALSE else pval
    }
    if(pval) {
      pvalue <- node$info$p.value
      plab <- ifelse(pvalue < 10^(-3L),
                     paste("p <", 10^(-3L)),
                     paste("p =", round(pvalue, digits = 3L)))
    } else {
      plab <- ""
    }
    return(c(varlab, plab))
  }
  
  maxstr <- function(node) {
    lab <- extract_label(node)
    klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
    lab <- c(lab, klab)
    lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
    lab <- lab[which.max(nchar(lab))]
    if(length(lab) < 1L) lab <- ""
    return(lab)
  }
  
  nstr <- maxstr(node_party(obj))
  if(nchar(nstr) < 6) nstr <- "aAAAAa"
  
  ### panel function for the inner nodes
  rval <- function(node, ...) {  
    
    # extract ...
    kwargs <- list(...)
    if(!is.null(kwargs$remove.nobs)) {remove.nobs <- kwargs$remove.nobs} else {remove.nobs <- FALSE}
    
    node_vp <- viewport(
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "strwidth", nstr)+unit(.2, "lines"), 
      height = unit(3, "lines"),
      name = paste("node_inner.left_rotated_tree", id_node(node), sep = ""),
      gp = gp
    )
    pushViewport(node_vp)
    
    xell <- c(seq(0, 0.2, by = 0.01),
              seq(0.2, 0.8, by = 0.05),
              seq(0.8, 1, by = 0.01))
    yell <- sqrt(xell * (1-xell))
    
    lab <- extract_label(node)
    
    #fill <- rep(fill, length.out = 2L)
    # grid.polygon(x = unit(c(xell, rev(xell)), "npc"),#x = unit(0.1, "npc"), y = unit(0.1, "npc"),
    #              y = unit(c(yell, -yell)+0.5, "npc"),
    #              gp = gpar(fill = fill[1], col=fill[1]))
    
    grid.rect(height=unit(0.15,"lines") ,gp = gpar(fill = fill, col=fill))

    ## FIXME: something more general instead of pval ?
    grid.text(lab[1L], y = unit(1.5 + 0.5 * (lab[2L] != ""), "lines")) #That's the variable name
    #grid.text('Test', y = unit(1, "lines")) #That's the variable name
    if(lab[2L] != "") grid.text(lab[2L], y = unit(1, "lines")) #Printing p-value
    
    
    if(id) {
      if(id_node(node)==1) {
        nodeIDvp <- viewport(x = unit(.5, "npc")+unit(.5, "strwidth", nstr)+unit(.1,"lines"), 
                             y = unit(.5, "npc"),
                             width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
                             height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])),
                             just="left")
        pushViewport(nodeIDvp)
        grid.rect(gp = gpar(fill = fill))
        grid.text(nam[id_node(node)]) #print node number
        popViewport()
        
        # Print number of observations
        
        if(!remove.nobs) {
          dat <- data_party(obj, id_node(node))
          yn <- dat[["(response)"]]
          wn <- dat[["(weights)"]]
          if(is.null(wn)) wn <- rep(1, NROW(yn))
          nodeNum <- viewport(x = unit(.5, "npc"),
                              #- unit(.5, "strwidth", nstr) + unit(0.2, "lines")+ max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])), 
                              y = unit(1, "npc"),
                              width = max(unit(.8, "lines"), unit(1.3, "strwidth", sprintf("n = %s", sum(wn)))),
                              height = max(unit(1, "lines"), unit(1.3, "strheight", sprintf("n = %s", sum(wn)))))
          pushViewport(nodeNum)
          # if (debug) {
          #   grid.rect(gp = gpar(fill = fill[2], lty="dotted", col="red"))
          # } else {
          grid.rect(gp = gpar(fill = fill[2], col="white"))#
          # }
          grid.text(sprintf("n = %s", sum(wn)), just = c("center", "center")) #print node number
          popViewport()
        }
      }
    }
    upViewport()
  }
  
  return(rval)
}
class(node_inner.left_rotated_tree) <- "grapcon_generator"

# ### Plot function for terminal node ###
# # To change it, write a new function and execute
# # R> plot.left_rotated_tree(..., terminal_panel = YOURFUNCTION)
# node_ecdf.left_rotated_tree <- function(obj, col = "black", bg = "white", ylines = 2,
#                                       id = TRUE, mainlab = NULL, gp = gpar())
# {
#   
#   ## extract response
#   y <- obj$fitted[["(response)"]]
#   stopifnot(inherits(y, "numeric") || inherits(y, "integer"))
#   
#   dostep <- function(f) {
#     x <- knots(f)
#     y <- f(x)
#     ### create a step function based on x, y coordinates
#     ### modified from `survival:print.survfit'
#     if (is.na(x[1] + y[1])) {
#       x <- x[-1]
#       y <- y[-1]
#     }
#     n <- length(x)
#     if (n > 2) {  
#       # replace verbose horizonal sequences like
#       # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
#       # with (1, .2), (3, .1).  They are slow, and can smear the looks
#       # of the line type.
#       dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
#       n2 <- sum(dupy)
#       
#       #create a step function
#       xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
#       yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
#       RET <- list(x = xrep, y = yrep)
#     } else {
#       if (n == 1) {
#         RET <- list(x = x, y = y)
#       } else {
#         RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
#       }
#     }
#     return(RET)
#   }
#   
#   ### panel function for ecdf in nodes
#   rval <- function(node, ...) {
#     
#     # extract ...
#     #pass
#     kwargs <- list(...)
#     if(!is.null(kwargs$remove.xaxis)) {remove.xaxis <- kwargs$remove.xaxis} else {remove.xaxis <- FALSE}
#     if(!is.null(kwargs$cex)) {cex <- kwargs$cex} else {cex <- 1}
#     
#     ## extract data
#     nid <- id_node(node)
#     dat <- data_party(obj, nid)
#     yn <- dat[["(response)"]]
#     wn <- dat[["(weights)"]]
#     if(is.null(wn)) wn <- rep(1, NROW(yn))
#     
#     #defining helper function
#     .pred_ecdf <- function(y, w) {
#       if (length(y) == 0) return(NA)
#       iw <- as.integer(round(w))
#       if (max(abs(w - iw)) < sqrt(.Machine$double.eps)) {
#         y <- rep(y, w)
#         return(ecdf(y))
#       } else {
#         stop("cannot compute empirical distribution function with non-integer weights")
#       }
#     }
#     
#     ## get ecdf in node
#     f <- .pred_ecdf(yn, wn)
#     a <- dostep(f)
#     
#     
#     ## set up plot
#     yscale <- c(0, 1)
#     xscale <- range(y)
#     a$x <- c(xscale[1], a$x[1], a$x, xscale[2])
#     a$x <- a$x - min(a$x)
#     a$x <- a$x / max(a$x)
#     a$y <- c(0, 0, a$y, 1)
#     
#     if(remove.xaxis) {
#       top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 2,
#                                               widths = unit(c(ylines, 1),
#                                                             c("lines", "null")),
#                                               heights = unit(c(1.2, 1), c("lines", "null"))),
#                          width = unit(1, "npc"),
#                          height = unit(1, "npc"),
#                          name = paste("node_ecdf", nid, sep = ""), gp = gp)
#     } else {
#       top_vp <- viewport(y=unit(.5, "npc")+unit(1,"lines"),
#                          layout = grid.layout(nrow = 2, ncol = 2,
#                                               widths = unit(c(ylines, 1),
#                                                             c("lines", "null")),
#                                               heights = unit(c(1.2, 1), c("lines", "null"))),
#                          width = unit(1, "npc"),
#                          height = unit(1, "npc")-unit(2, "lines"),
#                          name = paste("node_ecdf", nid, sep = ""), gp = gp)
#     }
#     
#     pushViewport(top_vp)
#     grid.rect(gp = gpar(fill = bg, col = 0))
#     
#     ## number of observations
#     top <- viewport(layout.pos.col=2, layout.pos.row=1)
#     pushViewport(top)
#     
#     n.text <- sprintf("n = %s", sum(wn))
#     grid.rect(x = unit(0, "npc"), 
#               y = unit(.5, "lines"),
#               height = unit(1, "lines"), 
#               width = unit(1, "strwidth", n.text)+unit(.2, "lines"), 
#               gp = gpar(fill = "transparent", col = 1), 
#               just= "left")
#     grid.text(n.text, x = unit(0, "npc")+unit(.1, "lines"), 
#                       y = unit(.5, "lines"),
#                       just="left")
#     popViewport()
#     
#     plot <- viewport(layout.pos.col=2, layout.pos.row=2,
#                      xscale=xscale, yscale=yscale,
#                      name = paste0("node_surv", nid, "plot"),
#                      clip = FALSE)
#     
#     pushViewport(plot)
#     if(!remove.xaxis) {grid.xaxis()}
#     grid.yaxis()
#     grid.rect(gp = gpar(fill = "transparent"))
#     grid.clip()
#     grid.lines(a$x, a$y, gp = gpar(col = col))
#     
#     # Dashed line for mean
#     x.mean <- a$x[a$y[-length(a$y)]<=0.5 & a$y[-1]>0.5]
#     grid.lines(x=c(x.mean, x.mean), y=c(0,1), gp = gpar(col=2, lty="dashed"))
#     
#     upViewport(2)
#   }
#   
#   return(rval)
# }
# class(node_ecdf.left_rotated_tree) <- "grapcon_generator"


### Plot function for the edge panel ###
edge_simple.left_rotated_tree <- function(obj, digits = 3, abbreviate = FALSE,
                                        justmin = Inf, just = c("alternate", "increasing", "decreasing", "equal"),
                                        fill = "white")
{
  
  meta <- obj$data
  
  justfun <- function(i, split) {
    myjust <- if(mean(nchar(split)) > justmin) {
      match.arg(just, c("alternate", "increasing", "decreasing", "equal"))
    } else {
      "equal"
    }
    k <- length(split)
    rval <- switch(myjust,
                   "equal" = rep.int(0, k),
                   "alternate" = rep(c(0.5, -0.5), length.out = k),
                   "increasing" = seq(from = -k/2, to =  k/2, by = 1),
                   "decreasing" = seq(from =  k/2, to = -k/2, by = -1)
    )
    unit(0.5, "npc") + unit(rval[i], "lines")
  }
  
  ### panel function for simple edge labelling
  function(node, i, debug=FALSE) {
    split <- character_split(split_node(node), meta, digits = digits)$levels
    y <- justfun(i, split)
    split <- split[i]
    # try() because the following won't work for split = "< 10 Euro", for example.
    if(any(grep(">", split) > 0) | any(grep("<", split) > 0)) {
      tr <- suppressWarnings(try(parse(text = paste("phantom(0)", split)), silent = TRUE))
      if(!inherits(tr, "try-error")) split <- tr
    }
    grid.rect(x=0, y = y, gp = gpar(fill = fill, col = fill), width = unit(1, "strwidth", split),just="left")
    grid.text(split, x=0, y = y, just = "left")
  }
}
class(edge_simple.left_rotated_tree) <- "grapcon_generator"
