##### make_dist_list
# input: gamlss.dist family object
# builds a list with all necessary functions/informations for distfitlist

# distfit(y, family = make_dist_list(BI, bd = 10))
# distfit(y, family = BI, bd = 10)

#### TO DO: 

# seperate make_dist_list (as a function outside of distfitlist)
# bd can be an input argument in make_dist_list ( => bd is stored in the local environment), but once the list is returned, bd should not appear any more
# such that when writing a list (for distributions without bd) the parameter bd should not appear at any time
# wrap another function around gamlss.dist family object functions  such that bd is fixed inside this wrapper function and then the dldm etc. functions are called
# input of wrapper function mustn't include bd (make_dist_list is already the wrapper function for all included functions)

# within distfitlist: if hdist is available in the list, the hessian should be calculated analytically using hdist
# while building the list: if hdist is availabe, but the input argument type.hessian is set to "numeric" -> remove hdist



make_dist_list <- function(family, bd = 1) {
  
  
  ## list of families which require an additional parameter bd (binomial denominator)
  # by default bd is set to 1 
  .distfit.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial", "ZIBI", "ZIBB", "ZABI", "ZABB") # binomial denominators
  # if(any(family$family%in%.distfit.bi.list)) 

  np <- sum(family$parameter == TRUE)
  
  ## families which have fixed parameters: LNO (log normal (Box-Cox)), NET
  # define fixed parameters locally within make_dist_list, such that the fixed parameters don't appear in the output list anymore
  # for LNO a value. nu.start is needed for the evaluation of mu.initial
  if(family$family[1] == "LNO") nu.start <- NULL
  
  if(FALSE) {
    if(family$nopar != np){
      if(family$family[1] == "LNO") nu <- nu.start <- 0
      if(family$family[1] == "NET") {
        nu <- 1.5
        tau <- 2
      }
      
      #general form: (problem: y required in .initial)
      #if(!family$parameters$mu) {eval(family$mu.initial)}
      #if(!family$parameters$sigma) {eval(family$sigma.initial)}
      #if(!family$parameters$nu) {eval(family$nu.initial)}
      #if(!family$parameters$tau) {eval(family$tau.initial)}
    }
  }
  
  
  ## notation:
  # par ... distribution parameters (mu, sigma, nu, tau)
  # eta ... coefficients of the linear predictor, here: intercept (g(mu)=eta[1], g(sigma)=eta[2], g(nu)=eta[3], g(tau)=eta[4])
  
  # if(np > 0L) m <- family$mu.linkinv(eta[1L])          # m ... mu           eta[1] ... g(mu)        g ... link function
  # if(np > 1L) s <- family$sigma.linkinv(eta[2L])       # s ... sigma        eta[2] ... g(sigma)     g ... link function
  # if(np > 2L) v <- family$nu.linkinv(eta[3L])          # v ... nu           eta[3] ... g(nu)        g ... link function
  # if(np > 3L) t <- family$tau.linkinv(eta[4L])         # t ... tau          eta[4] ... g(tau)       g ... link function
  
  
  
  ## Define all necessary functions depending on the number of parameters
  
  # get parameters of a function f, return vector with the indices of the necessary input parameters 
  getpar <- function(f){
    arguments <- names(formals(f))
    par.id <- c()
    if("mu" %in% arguments) par.id <- c(par.id, 1)
    if("sigma" %in% arguments) par.id <- c(par.id, 2)
    if("nu" %in% arguments) par.id <- c(par.id, 3)
    if("tau" %in% arguments) par.id <- c(par.id, 4)
    
    #if(family$nopar != np){
    #  if((family$family[1] == "LNO") par.id <- par.id[par.id != 3]
    #  if(family$family[1] == "NET") par.id <- par.id[par.id != c(3,4)]
    #}
    
    return(par.id)
  }
  
  
  # get derivative function as required in the list (input: gamlss.dist derivative function)
  getderivfun <- function(fun){
    arg <- names(formals(fun))
    par.id <- getpar(fun)
    if("y" %in% arg) {
      if("bd" %in% arg) {
        derivfun <- function(y, par) {
          input <- list()
          input$y <- y
          input <- c(input, par[par.id])
          input$bd <- bd
          val <- do.call(fun, input)
          if(length(val) == 1L) val <- rep(val, length(y))
          return(val)
        }
      } else {
        derivfun <- function(y, par) {
          input <- list()
          input$y <- y
          input <- c(input, par[par.id])
          val <- do.call(fun, input)
          if(length(val) == 1L) val <- rep(val, length(y))
          return(val)
        }    
      } 
    } else {
      if("bd" %in% arg) {
        derivfun <- function(y, par) {
          input <- list()
          input <- c(input, par[par.id])
          input$bd <- bd
          return(rep(do.call(fun, input), length(y)))
        }
      } else {
        derivfun <- function(y, par) {
          input <- list()
          input <- c(input, par[par.id])
          return(rep(do.call(fun, input), length(y)))
        }    
      }
    }
    return(derivfun)
  }
  
  
  ## define inner and outer derivative functions
  
  if(np > 0L){
    
    # define names
    parnames <- c("mu")
    linknames <- c(family$mu.link)
    if(family$mu.link == "identity") etanames <- c("mu") else etanames <- paste0(family$mu.link, "(mu)")
    
    
    # inner derivative functions (dmdeta, d2mdeta2)
    dmdeta <- function(eta) return(family$mu.dr(eta[1]))
    if(family$mu.link=="identity") d2mdeta2 <- function(eta) return(0)
    if(family$mu.link=="log")      d2mdeta2 <- function(eta) return(exp(eta[1]))
    if(family$mu.link=="logit")    d2mdeta2 <- function(eta) return(exp(eta[1]) * (exp(eta[1])-1) / ((1+exp(eta[1]))^3)) 
    
    # outer derivative functions (dldm, d2ldm2)
    dldm <- getderivfun(family$dldm)
    d2ldm2 <- getderivfun(family$d2ldm2)
    
    ## FIX ME: alternative version, pro: shorter, contra: if-conditions inside the function 
    if(FALSE){
      arg.dldm <- names(formals(family$dldm))
      par.id.dldm <- getpar(family$dldm)
      dldm <- function(y, par){
        input <- list()
        if("y" %in% arg.dldm) input$y <- y
        input <- c(input, par[par.id.dldm])
        if("bd" %in% arg.dldm) input$bd <- bd
        return(do.call(family$dldm, input))
      }
    }
    
  }
    
  
  if(np > 1L){
    
    # define names
    parnames <- c(parnames,"sigma")
    linknames <- c(linknames, family$sigma.link)
    if(family$sigma.link == "identity") etanames <- c(etanames,"sigma") else etanames <- c(etanames, paste0(family$sigma.link, "(sigma)"))
    
    # inner derivative functions (dddeta, d2ddeta2)     
    dddeta <- function(eta) return(family$sigma.dr(eta[2]))
    if(family$sigma.link=="identity") d2ddeta2 <- function(eta) return(0)
    if(family$sigma.link=="log")      d2ddeta2 <- function(eta) return(exp(eta[2]))
    if(family$sigma.link=="logit")    d2ddeta2 <- function(eta) return(exp(eta[2]) * (exp(eta[2])-1) / ((1+exp(eta[2]))^3)) 
    
    # outer derivative functions (dldd, d2ldd2, d2ldmdd)
    dldd <- getderivfun(family$dldd)
    d2ldd2 <- getderivfun(family$d2ldd2)
    d2ldmdd <- getderivfun(family$d2ldmdd)

  }
  
  
  if(np > 2L){
    
    # define names
    parnames <- c(parnames,"nu")
    linknames <- c(linknames, family$nu.link)
    if(family$nu.link == "identity") etanames <- c(etanames,"nu") else etanames <- c(etanames, paste0(family$nu.link, "(nu)"))
    
    # inner derivative functions (dvdeta, d2vdeta2)
    dvdeta <- function(eta) return(family$nu.dr(eta[3]))
    if(family$nu.link=="identity") d2vdeta2 <- function(eta) return(0)
    if(family$nu.link=="log")      d2vdeta2 <- function(eta) return(exp(eta[3]))
    if(family$nu.link=="logit")    d2vdeta2 <- function(eta) return(exp(eta[3]) * (exp(eta[3])-1) / ((1+exp(eta[3]))^3)) 
    
    # outer derivatives (dldv, d2ldv2, d2ldmdv, d2ldddv)
    dldv <- getderivfun(family$dldv)
    d2ldv2 <- getderivfun(family$d2ldv2)
    d2ldmdv <- getderivfun(family$d2ldmdv)
    d2ldddv <- getderivfun(family$d2ldddv)
  }
  
  
  if(np > 3L){
    
    # define names
    parnames <- c(parnames,"tau")
    linknames <- c(linknames, family$tau.link)
    if(family$tau.link == "identity") etanames <- c(etanames,"tau") else etanames <- c(etanames, paste0(family$tau.link, "(tau)"))
    
    # note: in this case/section no adaption for families of the list .distfit.bi.list since none of them includes the 4th parameter tau
    
    # inner derivatives (dtdeta, d2tdeta2)    
    dtdeta <- function(eta) return(family$tau.dr(eta[4]))
    if(family$tau.link=="identity")  d2tdeta2 <- function(eta) return(0)
    if(family$tau.link=="log")       d2tdeta2 <- function(eta) return(exp(eta[4]))
    if(family$tau.link=="logit")     d2tdeta2 <- function(eta) return(exp(eta[4]) * (exp(eta[4])-1) / ((1+exp(eta[4]))^3)) 
    
    # outer derivatives (dldt, d2ldt2, d2ldmdt, d2ldddt, d2ldvdt)
    dldt <- getderivfun(family$dldt)
    d2ldt2 <- getderivfun(family$d2ldt2)
    d2ldmdt <- getderivfun(family$d2ldmdt)
    d2ldddt <- getderivfun(family$d2ldddt)
    d2ldvdt <- getderivfun(family$d2ldvdt)
    
  }
  
  
  
  ## define complete derivative functions dpardeta, d2pardeta2, dldpar, d2ldpar2 according to the number of parameters
  
  if(np == 1L){
    
    # define function for the calculation of initial values
    ## FIX ME ## use weights?
    startfun <- function(y, weights = NULL) {
      if(!is.null(weights)) y <- rep(y, round(weights))
      mu <- NULL
      eval(family$mu.initial)
      starteta <- family$mu.linkfun(mean(mu))
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]))
      names(eta) <- etanames
      return(eta)
    }
  
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list of matrices:
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par)))
      
      # d2matrix is of size (1*ny x 1) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (1x1) is stored in d2list
      
      d2list <- list()
      ny <- length(y)
      if(is.Surv(y)) ny <- dim(y)[1]
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i),]
      }
      
      return(d2list)
    }
  }
  
  
  if(np == 2L){
    
    # define function for the calculation of initial values
    startfun <- function(y, weights = NULL) {
      if(!is.null(weights)) y <- rep(y, round(weights))
      mu <- sigma <- NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      starteta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)))
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]), family$mu.linkfun(par[2]))
      names(eta) <- etanames
      return(eta)
    }
    
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list of matrices:
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par), dldd(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par), d2ldmdd(y, par)),
                        cbind(d2ldmdd(y, par), d2ldd2(y, par)))
      
      # d2matrix is of size (2*ny x 2) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (2x2) is stored in d2list
      
      d2list <- list()
      ny <- length(y)
      if(is.Surv(y)) ny <- dim(y)[1]
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i),]
      }
      
      return(d2list)
    }
  }
  
  
  if(np == 3L){
    
    # define function for the calculation of initial values
    startfun <- function(y, weights = NULL) {
      if(!is.null(weights)) y <- rep(y, round(weights))
      mu <- sigma <- nu <-  NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      eval(family$nu.initial)
      starteta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)))
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]), family$mu.linkfun(par[2]), family$mu.linkfun(par[3]))
      names(eta) <- etanames
      return(eta)
    }
    
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta), dvdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta), d2vdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list:
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par), dldd(y, par), dldv(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par), d2ldmdd(y, par), d2ldmdv(y, par)),
                        cbind(d2ldmdd(y, par), d2ldd2(y, par), d2ldddv(y, par)),
                        cbind(d2ldmdv(y, par), d2ldddv(y, par), d2ldv2(y, par)))
      
      # d2matrix is of size (3*ny x 3) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (3x3) is stored in d2list
      
      d2list <- list()
      ny <- length(y)
      if(is.Surv(y)) ny <- dim(y)[1]
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i),]
      }
      
      return(d2list)
    }
  }
  
  
  if(np == 4L){
    
    # define function for the calculation of initial values
    startfun <- function(y, weights = NULL) {
      if(!is.null(weights)) y <- rep(y, round(weights))
      mu <- sigma <- nu <- tau <- NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      eval(family$nu.initial)
      eval(family$tau.initial)
      starteta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)), family$tau.linkfun(mean(tau)))
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]), family$mu.linkfun(par[2]), family$mu.linkfun(par[3]), family$mu.linkfun(par[4]))
      names(eta) <- etanames
      return(eta)
    }
    
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]), family$tau.linkinv(eta[4]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta), dvdeta(eta), dtdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta), d2vdeta2(eta), d2tdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list :
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par), dldd(y, par), dldv(y, par), dldt(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par), d2ldmdd(y, par), d2ldmdv(y, par), d2ldmdt(y, par)),
                        cbind(d2ldmdd(y, par), d2ldd2(y, par), d2ldddv(y, par), d2ldddt(y, par)),
                        cbind(d2ldmdv(y, par), d2ldddv(y, par), d2ldv2(y, par), d2ldvdt(y, par)),
                        cbind(d2ldmdt(y, par), d2ldddt(y, par), d2ldvdt(y, par), d2ldt2(y, par)))
      
      # d2matrix is of size (4*ny x 4) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (4x4) is stored in d2list
      
      d2list <- list()
      ny <- length(y)
      if(is.Surv(y)) ny <- dim(y)[1]
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i, 3*ny+i),]
      }
      
      return(d2list)
    }
  }
  
  
  
  
  ddist <- function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {
    par <- linkinv(eta)
    
    # fixed parameters do not need to be added here, because in the density function dLNO and dNET nu (and tau) are by default set to to c(0) or c(1.5, 2) respectively
    #if(family$nopar != np){
    #  if(family$family[1] == "LNO") par <- c(par, 0)
    #  if(family$family[1] == "NET") par <- c(par, 1.5, 2)
    #}
    
    input <- list()
    inputnames <- c("x", parnames, "log")
    input[[1]] <- y
    input <- c(input, par)                          # <- rep.int(par[i-1], length(y))   (FIX?)
    if(any(family$family%in%.distfit.bi.list)) {
      input <- c(input, bd)      # additional parameter bd (binomial denominator for families in .distfit.bi.list)
      inputnames <- c("x", parnames, "bd", "log")
    }
    input <- c(input, log)
    names(input) <- inputnames
    eval <- do.call(get(paste0("d", family$family[[1]])), input)
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      eval <- sum(weights * eval)
    }
    return(eval)
  }

  
  ## additional functions pdist, qdist, rdist
  if(FALSE) {
  if(any(family$family%in%.distfit.bi.list)){
    if(np == 1L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], bd = bd, log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], bd = bd, log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], bd = bd, log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1], bd = bd)
    }
    if(np == 2L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], sigma = par[2], bd = bd, log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], sigma = par[2], bd = bd, log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], sigma = par[2], bd = bd, log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1], sigma = par[2], bd = bd)
    }
    if(np == 3L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], sigma = par[2], nu = par[3], bd = bd, log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], sigma = par[2], nu = par[3], bd = bd, log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], sigma = par[2], nu = par[3], bd = bd, log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1], sigma = par[2], nu = par[3], bd = bd)
    }
    if(np == 4L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], bd = bd, log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], bd = bd, log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], bd = bd, log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], bd = bd)
    }
  } else {
    if(np == 1L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1])
    }
    if(np == 2L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], sigma = par[2], log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], sigma = par[2], log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], sigma = par[2], log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1], sigma = par[2])
    }
    if(np == 3L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], sigma = par[2], nu = par[3], log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], sigma = par[2], nu = par[3], log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], sigma = par[2], nu = par[3], log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1], sigma = par[2], nu = par[3])
    }
    if(np == 4L) { 
      #ddist <- function(x, par, log = FALSE) get(paste0("d",family$family[1]))(x, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], log = log)
      pdist <- function(q, par, log.p = FALSE) get(paste0("p",family$family[1]))(q, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], log.p = log.p)
      qdist <- function(p, par, log.p = FALSE) get(paste0("q",family$family[1]))(p, mu = par[1], sigma = par[2], nu = par[3], tau = par[4], log.p = log.p)
      rdist <- function(n, par) get(paste0("r",family$family[1]))(n, mu = par[1], sigma = par[2], nu = par[3], tau = par[4])
    }
  }
  }
    
   
  
  
  ## score / estfun (first-order partial derivatives of the (positive) log-likelihood function)
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {
    par <- linkinv(eta) 
    
    if(family$nopar != np){
      if(family$family[1] == "LNO") {
        par <- c(par, 0)
        eta <- c(eta, 0)
      }
      if(family$family[1] == "NET") {
        par <- c(par, 1.5, 2)
        eta <- c(eta, 1.5, 2)
      }
    }
    score <- t(t(dldpar(y, par)) * dpardeta(eta))
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  ## hessian (second-order partial derivatives of the (positive) log-likelihood function)
  hdist <- function(y, eta, weights = NULL) {    
    ny <- length(y)
    if(is.Surv(y)) ny <- dim(y)[1]
    if(is.null(weights)) weights <- rep.int(1, ny)
    
    par <- linkinv(eta)
    
    if(family$nopar != np){
      if(family$family[1] == "LNO") {
        par <- c(par, 0)
        eta <- c(eta, 0)
      }
      if(family$family[1] == "NET") {
        par <- c(par, 1.5, 2)
        eta <- c(eta, 1.5, 2)
      }
    }
    
    ## calculate derivative vectors / matrices / lists
    d2ldpar2.list <- d2ldpar2(y, par)
    dldpar.mat <- dldpar(y, par)
    
    dpardeta.vec <- dpardeta(eta)
    d2pardeta2.vec <- d2pardeta2(eta)
    
    ## calculation is split up in 2 parts: 
    # second outer derivatives times first inner derivatives and a diagonal matrix with the first outer and the second inner derivatives
    hess.list <- list()
    length(hess.list) <- length(d2ldpar2.list)
    for(i in 1:ny){
      hess.list[[i]] <- weights[i] * (t(d2ldpar2.list[[i]] * dpardeta.vec) * dpardeta.vec + diag(np) * as.vector(dldpar.mat[i,]) * d2pardeta2.vec)
    }
    
    ## calculate the sum over all matrices in the list (each for one observation)  
    hess <- Reduce('+', hess.list)
    hess <- as.matrix(hess)
    colnames(hess) <- rownames(hess) <-  etanames
    return(hess)
  }

  
  link <- linknames      # as defined above (within if(np == ))
  
  #linkfun <- linkfun     # as defined above (within if(np == ))
  
  #linkinv <- linkinv     # as defined above (within if(np == ))
  
  linkinvdr <- dpardeta

  mle <- FALSE

  
  dist_list <- list(family.name = paste(family$family[2], "Distribution", sep = " "),
                    ddist = ddist, 
                    sdist = sdist, 
                    hdist = hdist, 
                    link = link, 
                    linkfun = linkfun, 
                    linkinv = linkinv, 
                    linkinvdr = linkinvdr,
                    startfun = startfun,
                    mle = mle
                    )
}




###### dist_list for normal distribution
if(FALSE) {
  
  dist_list_normal <- list()
  
  parnames <- c("mu", "sigma")
  etanames <- c("mu", "log(sigma)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    par <- c(eta[1], exp(eta[2]))
    val <- 1/sqrt(2*pi*par[2]^2) * exp(- (y-par[1])^2 / (2*par[2]^2))
    if(log) val <- log(val)
    # val <- dnorm(y, mean = par[1], sd = par[2], log = log)
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      val <- sum(weights * val)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
    par <- c(eta[1], exp(eta[2]))                           
    score <- cbind(1/par[2]^2 * (y-par[1]), (-1/par[2] + ((y - par[1])^2)/(par[2]^3)) * exp(eta[2]))
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {    
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, ny)
    
    par <- c(eta[1], exp(eta[2]))                           
    
    d2ld.etamu2 <- sum(weights * rep.int(-1/par[2]^2, ny))
    d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-par[1])/par[2]^2)          # should be 0 for exact parameters (here ~ e-17 due to calculations)
    d2ld.etasigma2 <- sum(weights * (-2)*(y-par[1])^2/par[2]^2)         
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  

  startfun <- function(y, weights = NULL){
    if(!is.null(weights)) y <- rep(y, round(weights))
    ny <- length(y)
    mu <- mean(y)
    sigma <- sqrt(1/ny * sum((y - mu)^2))
    starteta <- c(mu, log(sigma))
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- TRUE
  
  dist_list_normal <- list(family.name = "Normal Distribution",
                           ddist = ddist, 
                           sdist = sdist, 
                           hdist = hdist, 
                           link = link, 
                           linkfun = linkfun, 
                           linkinv = linkinv, 
                           linkinvdr = linkinvdr,
                           startfun = startfun,
                           mle = mle
  )
}







###### dist_list for Poisson distribution
if(FALSE) {
  
  dist_list_poisson <- list()
  
  parnames <- c("mu")
  etanames <- c("log(mu)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    par <- exp(eta)
    val <- par^y / gamma(y+1) * exp(-par)
    if(log) val <- log(val)
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      val <- sum(weights * val)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
    par <- exp(eta)                           
    score <- (y - par)
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {
    if(is.null(weights)) weights <- rep.int(1, length(y))
    
    par <- exp(eta)                           
    hess <- rep(-par, length(y))
    hess <- as.matrix(sum(weights * hess))
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  link <- c("log")
  
  linkfun <- function(par) {
    eta <- log(par)
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- exp(eta)
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    if(!is.null(weights)) y <- rep(y, round(weights))
    ny <- length(y)
    mu <- mean(y)
    starteta <- log(mu)
    names(starteta) <- etanames
    return(starteta)
  }
  
  
  mle <- TRUE
  
  
  dist_list_poisson <- list(family.name = "Poisson Distribution",
                            ddist = ddist, 
                            sdist = sdist, 
                            hdist = hdist, 
                            link = link, 
                            linkfun = linkfun, 
                            linkinv = linkinv, 
                            linkinvdr = linkinvdr,
                            startfun = startfun,
                            mle = mle
  )
}






###### dist_list for exponential distribution
if(FALSE) {
  
  dist_list_exp <- list()
  
  parnames <- c("lambda")
  etanames <- c("log(lambda)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    par <- exp(eta)
    val <- par * exp(-par * y)
    if(log) val <- log(val)
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      val <- sum(weights * val)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {
    par <- exp(eta)                           
    score <- 1 - y * par
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights)) weights <- rep.int(1, length(y))
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {    
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, length(y))
    
    par <- exp(eta)                           
    hess <- -y * par
    hess <- as.matrix(sum(weights * hess))
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  link <- c("log")
  
  linkfun <- function(par) {
    eta <- log(par)
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- exp(eta)
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    if(!is.null(weights)) y <- rep(y, round(weights))
    ny <- length(y)
    lambda <- ny / sum(y)
    starteta <- log(lambda)
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- TRUE
  
  
  dist_list_exp <- list(family.name = "Exponential Distribution",
                        ddist = ddist, 
                        sdist = sdist, 
                        hdist = hdist, 
                        link = link, 
                        linkfun = linkfun, 
                        linkinv = linkinv, 
                        linkinvdr = linkinvdr,
                        startfun = startfun,
                        mle = mle
  )
}
