.prepare_args <- function(object, data, zformula, control, ...) {
  
  if (is.null(modcall <- getCall(object))) 
    stop("Need an object with call component, see getCall")
  
  ## get arguments for cforest call
  args <- list(...)
  # args$ntree <- ntree
  args$control <- control
  
  # ## arguments used in model
  # modargs <- as.list(modcall)[-1]
  
  ## formula and data
  if(is.null(data)) data <- eval(modcall$data)
  args$data <- data
  modformula <- eval(modcall$formula)
  
  ## in case I switch to mob
  # if(is.null(zformula)) zformula <- formula(~ .)
  # mobformula <- as.Formula(modformula, zformula)
  
  ## cforest formula
  if(is.null(zformula)) zformula <- "~ ."
  if(!is.character(zformula)) zformula <- paste(as.character(zformula), collapse = " ")
  modvars <- all.vars(modformula)
  args$formula <- as.formula(
    paste(paste(modvars, collapse = " + "), zformula)
  )
  
  return(args)
}