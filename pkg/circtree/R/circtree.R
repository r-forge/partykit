# -------------------------------------------------------------------
# - NAME:   circtree.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-07-31
# -------------------------------------------------------------------
# - PURPOSE: Wrapper function for distextree 
# -------------------------------------------------------------------

circtree <- function(formula,
                       data,
                       subset,
                       na.action = na.pass,
                       weights,
                       offset,
                       cluster,
                       control = distextree_control(...),
                       converged = NULL,
                       scores = NULL,
                       doFit = TRUE,
                       ...) {
  cl <- match.call()
  cl2 <- cl
  cl2[[1]] <- quote(disttree::distextree)
  cl2$family <- dist_vonmises()

  tree <- eval(cl2)
  tree$info$call <- cl

  class(tree) <- c("circtree", class(tree))
  tree
}

