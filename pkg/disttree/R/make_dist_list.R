##### make_dist_list
# input: gamlss.dist family object
# builds a list with all necessary functions/informations for distfit

make_dist_list <- function(family) {
  
  
  ## list of families which require an additional parameter bd (binomial denominator)
  # by default bd is set to 1
  .distfit.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial", "ZIBI", "ZIBB", "ZABI", "ZABB") # binomial denominators
  # if(any(family$family%in%.distfit.bi.list)) 

  np <- sum(family$parameter == TRUE)
  
  
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
    
    ## remove "bd" from the list of parameters in case it is included
    arg.bd <- FALSE
    if(any(family$family%in%.distfit.bi.list)){
      if("bd" %in% arguments){ 
        arg.bd <- TRUE
        arguments <- arguments[-(match("bd", arguments))]
      }
    }
    
    
    ## 2 cases: with or without y as input
    
    ## with y
    # if y is one of the input arguments, the derivative function returns a vector of the length of y
    
    if(arguments[1]=="y"){
      # 4 Parameter
      if(length(arguments)==5L) par.id <- c(0,1,2,3,4)        # f(y, mu=m, sigma=s, nu=v, tau=t)
      
      # 3 Parameter
      if(length(arguments)==4L){
        if(arguments[2]=="mu"){
          if(arguments[3]=="sigma"){
            if(arguments[4]=="nu"){
              par.id <- c(0,1,2,3)                            # f(y, mu=m, sigma=s, nu=v)
            } else {
              par.id <- c(0,1,2,4)                            # f(y, mu=m, sigma=s, tau=t)
            }
          } else {
            par.id <- c(0,1,3,4)                              # f(y, mu=m, nu=v, tau=t)
          }
        } else {
          par.id <- c(0,2,3,4)                                # f(y, sigma=s, nu=v, tau=t)
        }
      }
      
      # 2 Parameter
      if(length(arguments)==3L){
        if(arguments[2]=="mu"){
          if(arguments[3]=="sigma") par.id <- c(0,1,2)        # f(y, mu=m, sigma=s)
          if(arguments[3]=="nu")    par.id <- c(0,1,3)        # f(y, mu=m, nu=v)
          if(arguments[3]=="tau")   par.id <- c(0,1,4)        # f(y, mu=m, tau=t)
        }
        if(arguments[2]=="sigma"){
          if(arguments[3]=="nu")    par.id <- c(0,2,3)        # f(y, sigma=s, nu=v)
          if(arguments[3]=="tau")   par.id <- c(0,2,4)        # f(y, sigma=s, tau=t)
        }
        if(arguments[2]=="nu")      par.id <- c(0,3,4)        # f(y, nu=v, tau=t)
      }
      
      # 1 Parameter
      if(length(arguments)==2L){
        if(arguments[2]=="mu")      par.id <- c(0,1)          # f(y, mu=m)
        if(arguments[2]=="sigma")   par.id <- c(0,2)          # f(y, sigma=s)
        if(arguments[2]=="nu")      par.id <- c(0,3)          # f(y, nu=v)
        if(arguments[2]=="tau")     par.id <- c(0,4)          # f(y, tau=t)
      }
      
      # 0 Parameter
      if(length(arguments)==1L)     par.id <- c(0)            # f(y)
    } else {
      
      ## without y
      # in this case the functions return only a single value -> create vector by replicating this value ny times (necessary for the matrix (using sum and *1/ny)) 
      
      # 4 Parameter
      if(length(arguments)==4L) par.id <- c(1,2,3,4)          # f(mu=m, sigma=s, nu=v, tau=t)
      
      # 3 Parameter
      if(length(arguments)==3L){
        if(arguments[1]=="mu"){
          if(arguments[2]=="sigma"){
            if(arguments[3]=="nu"){
              par.id <- c(1,2,3)                              # f(mu=m, sigma=s, nu=v)
            } else {
              par.id <- c(1,2,4)                              # f(mu=m, sigma=s, tau=t)
            }
          } else {
            par.id <- c(1,3,4)                                # f(mu=m, nu=v, tau=t)
          }
        } else {
          par.id <- c(2,3,4)                                  # f(sigma=s, nu=v, tau=t)
        }
      }
      
      # 2 Parameter
      if(length(arguments)==2L){
        if(arguments[1]=="mu"){
          if(arguments[2]=="sigma") par.id <- c(1,2)          # f(mu=m, sigma=s)
          if(arguments[2]=="nu")    par.id <- c(1,3)          # f(mu=m, nu=v)
          if(arguments[2]=="tau")   par.id <- c(1,4)          # f(mu=m, tau=t)
        }
        if(arguments[1]=="sigma"){ 
          if(arguments[2]=="nu")    par.id <- c(2,3)          # f(sigma=s, nu=v)
          if(arguments[2]=="tau")   par.id <- c(2,4)          # f(sigma=s, tau=t)
        }
        if(arguments[1]=="nu")      par.id <- c(3,4)          # f(nu=v, tau=t)
      }
      
      # 1 Parameter
      if(length(arguments)==1L){
        if(arguments[1]=="mu")    par.id <- c(1)              # f(mu=m)
        if(arguments[1]=="sigma") par.id <- c(2)              # f(sigma=s)
        if(arguments[1]=="nu")    par.id <- c(3)              # f(nu=v)
        if(arguments[1]=="tau")   par.id <- c(4)              # f(tau=t)
      }
      
      # 0 Parameter
      if(length(arguments)==0L)  par.id <- NULL       ## fix: possible case? of which class is f when the derivative is a constant?
    }
    
    ## attach index for the binomial denominator parameter bd if it has been removed earlier
    if(arg.bd) par.id <- c(par.id,5)
    
    return(par.id)
  }
  
  
  
  ## define inner and outer derivative functions
  
  if(np > 0L){
    
    # define names
    par.names <- c("mu")
    link.names <- c(family$mu.link)
    if(family$mu.link == "identity") eta.names <- c("mu") else eta.names <- paste0(family$mu.link, "(mu)")
    
    # inner derivative functions (dmdeta, d2mdeta2)
    dmdeta <- function(eta) return(family$mu.dr(eta[1]))
    if(family$mu.link=="identity") d2mdeta2 <- function(eta) return(0)
    if(family$mu.link=="log")      d2mdeta2 <- function(eta) return(exp(eta[1]))
    if(family$mu.link=="logit")    d2mdeta2 <- function(eta) return(exp(eta[1]) * (exp(eta[1])-1) / ((1+exp(eta[1]))^3)) 
    
    # outer derivative functions (dldm, d2ldm2)
    par.id.dldm <- getpar(family$dldm)
    if(par.id.dldm[1] == 0L) {
      if(length(par.id.dldm) == 1L){
        dldm <- function(y, par) return(family$dldm(y))
      } else {
        dldm <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.dldm){
            for (i in 2:(length(par.id.dldm)-1)) input[[i]] <- rep.int(par[par.id.dldm[i]], ny)
            input[[length(par.id.dldm)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.dldm)) input[[i]] <- rep.int(par[par.id.dldm[i]], ny)
          }
          return(do.call(family$dldm, input))
        }
      }
    } else {
      dldm <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.dldm){
          for (i in 1:(length(par.id.dldm)-1)) input[[i]] <- par[par.id.dldm[i]]
          input[[length(par.id.dldm)]] <- bd
        } else {
          for (i in 1:length(par.id.dldm)) input[[i]] <- par[par.id.dldm[i]]
        }
        return(rep.int(do.call(family$dldm, input), ny))
      }
    }
    
    par.id.d2ldm2 <- getpar(family$d2ldm2)
    if(par.id.d2ldm2[1] == 0L){
      if(length(par.id.d2ldm2) == 1L){
        d2ldm2 <- function(y, par) return(family$d2ldm2(y))
      } else {
        d2ldm2 <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.d2ldm2){
            for (i in 2:(length(par.id.d2ldm2)-1)) input[[i]] <- rep.int(par[par.id.d2ldm2[i]], ny)
            input[[length(par.id.d2ldm2)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.d2ldm2)) input[[i]] <- rep.int(par[par.id.d2ldm2[i]], ny)
          }
          return(do.call(family$d2ldm2, input))
        }
      }
    } else {
      d2ldm2 <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.d2ldm2){
          for (i in 1:(length(par.id.d2ldm2)-1)) input[[i]] <- par[par.id.d2ldm2[i]]
          input[[length(par.id.d2ldm2)]] <- bd
        } else {
          for (i in 1:length(par.id.d2ldm2)) input[[i]] <- par[par.id.d2ldm2[i]]
        }
        return(rep.int(do.call(family$d2ldm2, input), ny))
      }
    }
  }
  
  
  if(np > 1L){
    
    # define names
    par.names <- c(par.names,"sigma")
    link.names <- c(link.names, family$sigma.link)
    if(family$sigma.link == "identity") eta.names <- c(eta.names,"sigma") else eta.names <- c(eta.names, paste0(family$sigma.link, "(sigma)"))
    
    # inner derivative functions (dddeta, d2ddeta2)     
    dddeta <- function(eta) return(family$sigma.dr(eta[2]))
    if(family$sigma.link=="identity") d2ddeta2 <- function(eta) return(0)
    if(family$sigma.link=="log")      d2ddeta2 <- function(eta) return(exp(eta[2]))
    if(family$sigma.link=="logit")    d2ddeta2 <- function(eta) return(exp(eta[2]) * (exp(eta[2])-1) / ((1+exp(eta[2]))^3)) 
    
    # outer derivative functions (dldd, d2ldd2, d2ldmdd)
    par.id.dldd <- getpar(family$dldd)
    if(par.id.dldd[1] == 0L){
      if(length(par.id.dldd) == 1L){
        dldd <- function(y, par) return(family$dldd(y))
      } else {
        dldd <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.dldd){
            for (i in 2:(length(par.id.dldd)-1)) input[[i]] <- rep.int(par[par.id.dldd[i]], ny)
            input[[length(par.id.dldd)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.dldd)) input[[i]] <- rep.int(par[par.id.dldd[i]], ny)
          }
          return(do.call(family$dldd, input))
        }
      }
    } else {
      dldd <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.dldd){
          for (i in 1:(length(par.id.dldd)-1)) input[[i]] <- par[par.id.dldd[i]]
          input[[length(par.id.dldd)]] <- bd
        } else {
          for (i in 1:length(par.id.dldd)) input[[i]] <- par[par.id.dldd[i]]
        }
        return(rep.int(do.call(family$dldd, input), ny))
      }
    }
    
    par.id.d2ldd2 <- getpar(family$d2ldd2)
    if(par.id.d2ldd2[1] == 0L){
      if(length(par.id.d2ldd2) == 1L){
        d2ldd2 <- function(y, par) return(family$d2ldd2(y))
      } else {
        d2ldd2 <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.d2ldd2){
            for (i in 2:(length(par.id.d2ldd2)-1)) input[[i]] <- rep.int(par[par.id.d2ldd2[i]], ny)
            input[[length(par.id.d2ldd2)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.d2ldd2)) input[[i]] <- rep.int(par[par.id.d2ldd2[i]], ny)
          }
          return(do.call(family$d2ldd2, input))
        }
      }
    } else {
      d2ldd2 <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.d2ldd2){
          for (i in 1:(length(par.id.d2ldd2)-1)) input[[i]] <- par[par.id.d2ldd2[i]]
          input[[length(par.id.d2ldd2)]] <- bd
        } else {
          for (i in 1:length(par.id.d2ldd2)) input[[i]] <- par[par.id.d2ldd2[i]]
        }
        return(rep.int(do.call(family$d2ldd2, input), ny))
      }
    }
    
    par.id.d2ldmdd <- getpar(family$d2ldmdd)
    if(par.id.d2ldmdd[1] == 0L){
      if(length(par.id.d2ldmdd) == 1L){
        d2ldmdd <- function(y, par) return(family$d2ldmdd(y))
      } else {
        d2ldmdd <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.d2ldmdd){
            for (i in 2:(length(par.id.d2ldmdd)-1)) input[[i]] <- rep.int(par[par.id.d2ldmdd[i]], ny)
            input[[length(par.id.d2ldmdd)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.d2ldmdd)) input[[i]] <- rep.int(par[par.id.d2ldmdd[i]], ny)
          }
          return(do.call(family$d2ldmdd, input))
        }
      }
    } else {
      d2ldmdd <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.d2ldmdd){
          for (i in 1:(length(par.id.d2ldmdd)-1)) input[[i]] <- par[par.id.d2ldmdd[i]]
          input[[length(par.id.d2ldmdd)]] <- bd
        } else {
          for (i in 1:length(par.id.d2ldmdd)) input[[i]] <- par[par.id.d2ldmdd[i]]
        }
        return(rep.int(do.call(family$d2ldmdd, input), ny))
      }
    }
  }
  
  
  if(np > 2L){
    
    # define names
    par.names <- c(par.names,"nu")
    link.names <- c(link.names, family$nu.link)
    if(family$nu.link == "identity") eta.names <- c(eta.names,"nu") else eta.names <- c(eta.names, paste0(family$nu.link, "(nu)"))
    
    # inner derivative functions (dvdeta, d2vdeta2)
    dvdeta <- function(eta) return(family$nu.dr(eta[3]))
    if(family$nu.link=="identity") d2vdeta2 <- function(eta) return(0)
    if(family$nu.link=="log")      d2vdeta2 <- function(eta) return(exp(eta[3]))
    if(family$nu.link=="logit")    d2vdeta2 <- function(eta) return(exp(eta[3]) * (exp(eta[3])-1) / ((1+exp(eta[3]))^3)) 
    
    # outer derivatives (dldv, d2ldv2, d2ldmdv, d2ldddv)
    par.id.dldv <- getpar(family$dldv)
    if(par.id.dldv[1] == 0L){
      if(length(par.id.dldv) == 1L){
        dldv <- function(y, par) return(family$dldv(y))
      } else {
        dldv <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.dldv){
            for (i in 2:(length(par.id.dldv)-1)) input[[i]] <- rep.int(par[par.id.dldv[i]], ny)
            input[[length(par.id.dldv)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.dldv)) input[[i]] <- rep.int(par[par.id.dldv[i]], ny)
          }
          return(do.call(family$dldv, input))
        }
      }
    } else {
      dldv <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.dldv){
          for (i in 1:(length(par.id.dldv)-1)) input[[i]] <- par[par.id.dldv[i]]
          input[[length(par.id.dldv)]] <- bd
        } else {
          for (i in 1:length(par.id.dldv)) input[[i]] <- par[par.id.dldv[i]]
        }
        return(rep.int(do.call(family$dldv, input), ny))
      }
    }
    
    par.id.d2ldv2 <- getpar(family$d2ldv2)
    if(par.id.d2ldv2[1] == 0L){
      if(length(par.id.d2ldv2) == 1L){
        d2ldv2 <- function(y, par) return(family$d2ldv2(y))
      } else {
        d2ldv2 <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.d2ldv2){
            for (i in 2:(length(par.id.d2ldv2)-1)) input[[i]] <- rep.int(par[par.id.d2ldv2[i]], ny)
            input[[length(par.id.d2ldv2)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.d2ldv2)) input[[i]] <- rep.int(par[par.id.d2ldv2[i]], ny)
          }
          return(do.call(family$d2ldv2, input))
        }
      }
    } else {
      d2ldv2 <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.d2ldv2){
          for (i in 1:(length(par.id.d2ldv2)-1)) input[[i]] <- par[par.id.d2ldv2[i]]
          input[[length(par.id.d2ldv2)]] <- bd
        } else {
          for (i in 1:length(par.id.d2ldv2)) input[[i]] <- par[par.id.d2ldv2[i]]
        }
        return(rep.int(do.call(family$d2ldv2, input), ny))
      }
    }
    
    par.id.d2ldmdv <- getpar(family$d2ldmdv)
    if(par.id.d2ldmdv[1] == 0L) { 
      if(length(par.id.d2ldmdv) == 1L){
        d2ldmdv <- function(y, par) return(family$d2ldmdv(y))
      } else {
        d2ldmdv <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.d2ldmdv){
            for (i in 2:(length(par.id.d2ldmdv)-1)) input[[i]] <- rep.int(par[par.id.d2ldmdv[i]], ny)
            input[[length(par.id.d2ldmdv)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.d2ldmdv)) input[[i]] <- rep.int(par[par.id.d2ldmdv[i]], ny)
          }
          return(do.call(family$d2ldmdv, input))
        }
      }
    } else {
      d2ldmdv <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.d2ldmdv){
          for (i in 1:(length(par.id.d2ldmdv)-1)) input[[i]] <- par[par.id.d2ldmdv[i]]
          input[[length(par.id.d2ldmdv)]] <- bd
        } else {
          for (i in 1:length(par.id.d2ldmdv)) input[[i]] <- par[par.id.d2ldmdv[i]]
        }
        return(rep.int(do.call(family$d2ldmdv, input), ny))
      }
    }
    
    par.id.d2ldddv <- getpar(family$d2ldddv)
    if(par.id.d2ldddv[1] == 0L) { 
      if(length(par.id.d2ldddv) == 1L){
        d2ldddv <- function(y, par) return(family$d2ldddv(y))
      } else {
        d2ldddv <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          if(5%in%par.id.d2ldddv){
            for (i in 2:(length(par.id.d2ldddv)-1)) input[[i]] <- rep.int(par[par.id.d2ldddv[i]], ny)
            input[[length(par.id.d2ldddv)]] <- rep.int(bd, ny)
          } else {
            for (i in 2:length(par.id.d2ldddv)) input[[i]] <- rep.int(par[par.id.d2ldddv[i]], ny)
          }
          return(do.call(family$d2ldddv, input))
        }
      }
    } else {
      d2ldddv <- function(y, par) {
        ny <- length(y)
        input <- list()
        if(5%in%par.id.d2ldddv){
          for (i in 1:(length(par.id.d2ldddv)-1)) input[[i]] <- par[par.id.d2ldddv[i]]
          input[[length(par.id.d2ldddv)]] <- bd
        } else {
          for (i in 1:length(par.id.d2ldddv)) input[[i]] <- par[par.id.d2ldddv[i]]
        }
        return(rep.int(do.call(family$d2ldddv, input), ny))
      }
    }
  }
  
  
  if(np > 3L){
    
    # define names
    par.names <- c(par.names,"tau")
    link.names <- c(link.names, family$tau.link)
    if(family$tau.link == "identity") eta.names <- c(eta.names,"tau") else eta.names <- c(eta.names, paste0(family$tau.link, "(tau)"))
    
    # note: in this case/section no adaption for families of the list .distfit.bi.list since none of them includes the 4th parameter tau
    
    # inner derivatives (dtdeta, d2tdeta2)    
    dtdeta <- function(eta) return(family$tau.dr(eta[4]))
    if(family$tau.link=="identity")  d2tdeta2 <- function(eta) return(0)
    if(family$tau.link=="log")       d2tdeta2 <- function(eta) return(exp(eta[4]))
    if(family$tau.link=="logit")     d2tdeta2 <- function(eta) return(exp(eta[4]) * (exp(eta[4])-1) / ((1+exp(eta[4]))^3)) 
    
    # outer derivatives (dldt, d2ldt2, d2ldmdt, d2ldddt, d2ldvdt)
    par.id.dldt <- getpar(family$dldt)
    if(par.id.dldt[1] == 0L){
      if(length(par.id.dldt) == 1L){
        dldt <- function(y, par) return(family$dldt(y))
      } else {
        dldt <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.dldt)) input[[i]] <- rep.int(par[par.id.dldt[i]], ny)
          return(do.call(family$dldt, input))
        }
      }
    } else {
      dldt <- function(y, par) {
        ny <- length(y)
        input <- list()
        for (i in 1:length(par.id.dldt)) input[[i]] <- par[par.id.dldt[i]]
        return(rep.int(do.call(family$dldt, input), ny))
      }
    }
    
    par.id.d2ldt2 <- getpar(family$d2ldt2)
    if(par.id.d2ldt2[1] == 0L){
      if(length(par.id.d2ldt2) == 1L){
        d2ldt2 <- function(y, par) return(family$d2ldt2(y))
      } else {
        d2ldt2 <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldt2)) input[[i]] <- rep.int(par[par.id.d2ldt2[i]], ny)
          return(do.call(family$d2ldt2, input))
        }
      }
    } else {
      d2ldt2 <- function(y, par) {
        ny <- length(y)
        input <- list()
        for (i in 1:length(par.id.d2ldt2)) input[[i]] <- par[par.id.d2ldt2[i]]
        return(rep.int(do.call(family$d2ldt2, input), ny))
      }
    }
    
    par.id.d2ldmdt <- getpar(family$d2ldmdt)
    if(par.id.d2ldmdt[1] == 0L) { 
      if(length(par.id.d2ldmdt) == 1L){
        d2ldmdt <- function(y, par) return(family$d2ldmdt(y))
      } else {
        d2ldmdt <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldmdt)) input[[i]] <- rep.int(par[par.id.d2ldmdt[i]], ny)
          return(do.call(family$d2ldmdt, input))
        }
      }
    } else {
      d2ldmdt <- function(y, par) {
        ny <- length(y)
        input <- list()
        for (i in 1:length(par.id.d2ldmdt)) input[[i]] <- par[par.id.d2ldmdt[i]]
        return(rep.int(do.call(family$d2ldmdt, input), ny))
      }
    } 
    
    par.id.d2ldddt <- getpar(family$d2ldddt)
    if(par.id.d2ldddt[1] == 0L) { 
      if(length(par.id.d2ldddt) == 1L){
        d2ldddt <- function(y, par) return(family$d2ldddt(y))
      } else {
        d2ldddt <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldddt)) input[[i]] <- rep.int(par[par.id.d2ldddt[i]], ny)
          return(do.call(family$d2ldddt, input))
        }
      }
    } else {
      d2ldddt <- function(y, par) {
        ny <- length(y)
        input <- list()
        for (i in 1:length(par.id.d2ldddt)) input[[i]] <- par[par.id.d2ldddt[i]]
        return(rep.int(do.call(family$d2ldddt, input), ny))
      }
    }
    
    par.id.d2ldvdt <- getpar(family$d2ldvdt)
    if(par.id.d2ldvdt[1] == 0L) { 
      if(length(par.id.d2ldvdt) == 1L){
        d2ldvdt <- function(y, par) return(family$d2ldvdt(y))
      } else {
        d2ldvdt <- function(y, par) {
          ny <- length(y)
          input <- list()
          input[[1]] <- y
          for (i in 2:length(par.id.d2ldvdt)) input[[i]] <- rep.int(par[par.id.d2ldvdt[i]], ny)
          return(do.call(family$d2ldvdt, input))
        }
      }
    } else {
      d2ldvdt <- function(y, par) {
        ny <- length(y)
        input <- list()
        for (i in 1:length(par.id.d2ldvdt)) input[[i]] <- par[par.id.d2ldvdt[i]]
        return(rep.int(do.call(family$d2ldvdt, input), ny))
      }
    }
  }
  
  
  
  ## define complete derivative functions dpardeta, d2pardeta2, dldpar, d2ldpar2 according to the number of parameters
  
  if(np == 1L){
    
    # define function for the calculation of initial values
    ## FIXME ## use weights?
    initialize <- function(y) {
      mu <- NULL
      eval(family$mu.initial)
      start.eta <- family$mu.linkfun(mean(mu))
      names(start.eta) <- eta.names
      return(start.eta)
    }
    
    
    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]))
      names(par) <- par.names
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
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i),]
      }
      
      return(d2list)
    }
  }
  
  
  if(np == 2L){
    
    # define function for the calculation of initial values
    initialize <- function(y) {
      mu <- sigma <- NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      start.eta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)))
      names(start.eta) <- eta.names
      return(start.eta)
    }
    
    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]))
      names(par) <- par.names
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
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i),]
      }
      
      return(d2list)
    }
  }
  
  
  if(np == 3L){
    
    # define function for the calculation of initial values
    initialize <- function(y) {
      mu <- sigma <- nu <-  NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      eval(family$nu.initial)
      start.eta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)))
      names(start.eta) <- eta.names
      return(start.eta)
    }
    
    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]))
      names(par) <- par.names
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
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i),]
      }
      
      return(d2list)
    }
  }
  
  
  if(np == 4L){
    
    # define function for the calculation of initial values
    initialize <- function(y) {
      mu <- sigma <- nu <- tau <- NULL
      eval(family$mu.initial)
      eval(family$sigma.initial)
      eval(family$nu.initial)
      eval(family$tau.initial)
      start.eta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)), family$tau.linkfun(mean(tau)))
      names(start.eta) <- eta.names
      return(start.eta)
    }
    
    # define function to get distribution parameters
    distpar <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]), family$tau.linkinv(eta[4]))
      names(par) <- par.names
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
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i, 3*ny+i),]
      }
      
      return(d2list)
    }
  }
  
  
  
  
  ddist <- function(y, par, log = TRUE, type = "parameter") {
    if(type == "link") par <- distpar(par)
    input <- list()
    input.names <- c("x", par.names, "log")
    input[[1]] <- y
    for(i in 2:(length(par)+1)) input[[i]] <- par[i-1]                           # <- rep.int(par[i-1], ny)   (FIX?)
    if(any(family$family%in%.distfit.bi.list)) {
      input[[length(par)+2]] <- bd      # additional parameter bd (binomial denominator for families in .distfit.bi.list)
      input.names <- c("x", par.names, "bd", "log")
    }
    input <- append(input, log)
    names(input) <- input.names
    eval <- do.call(get(paste0("d", family$family[[1]])), input)

    return(eval)
  }
  

  
  ## additional functions pdist, qdist, rdist
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
    
   
  
  
  ## score / estfun (for positive loglikelihood)
  sdist <- function(y, par, type = "parameter") {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    if(type == "parameter") {
      score <- dldpar(y, par)
      score <- as.matrix(score)
      colnames(score) <- par.names
    }
     if(type == "link") {
       eta <- par
       par <- distpar(eta)                           
       score <- t(t(dldpar(y, par)) * dpardeta(eta))
       score <- as.matrix(score)
       colnames(score) <- eta.names
     }
    return(score)
  }
  
  
  ## FIX ME: would it be better to return the list of hessian matrices (one matrix for each observation) ?
  ## hessian (for positive loglikelihood)
  hdist <- function(y, par, type = "parameter", weights = NULL) {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, ny)
    
    if(type == "parameter") {
      hess.list <- d2ldpar2(y, par)
      for(i in 1:ny){
        hess.list[[i]] <- weights[i] * hess.list[[i]]
      }
      hess <- Reduce('+', hess.list)
      # hess sum(weights * hess.list)
      hess <- as.matrix(hess)
      colnames(hess) <- rownames(hess) <-  par.names
    }
    
    if(type == "link") {
      eta <- par
      par <- distpar(eta)                           
      
      ## calculate derivative vectors / matrices / lists
      d2ldpar2.list <- d2ldpar2(y, par)
      dldpar.mat <- dldpar(y, par)
      
      dpardeta.vec <- dpardeta(eta)
      d2pardeta2.vec <- d2pardeta2(eta)
      
      ## calculation is split up in 2 parts: 
      # 2nd outer derivatives times first inner derivatives and a diagonal matrix with the first outer and the second inner derivatives
      hess.list <- list()
      length(hess.list) <- length(d2ldpar2.list)
      for(i in 1:ny){
        hess.list[[i]] <- weights[i] * (t(d2ldpar2.list[[i]] * dpardeta.vec) * dpardeta.vec + diag(np) * as.vector(dldpar.mat[i,]) * d2pardeta2.vec)
      }
      
      ## calculate the sum over all matrices in the list (each for one observation)  
      hess <- Reduce('+', hess.list)
      hess <- as.matrix(hess)
      colnames(hess) <- rownames(hess) <-  eta.names
    }
    
    return(hess)
    # return(hess, hess.list)
  }
  
  
  link <- link.names
  
  
  link.fun <- function(par) {
    eta <- NULL
    if(np > 0) eta[1] <- family$mu.linkfun(par[1])
    if(np > 1) eta[2] <- family$sigma.linkfun(par[2])
    if(np > 2) eta[3] <- family$nu.linkfun(par[3])
    if(np > 3) eta[4] <- family$tau.linkfun(par[4])
    names(eta) <- eta.names
    return(eta)
  }
  
  
  link.inv <- distpar
  
  
  link.inv.der <- dpardeta
  
  
  start.eta <- initialize
  
  
  start.par <- function(y) return(distpar(initialize(y)))
  

  use.optim <- TRUE
  
  if(use.optim) closed.mle <- NULL
  
  dist_list <- list(family.name = paste(family$family[2], "Distribution", sep = " "),
                    ddist = ddist, 
                    sdist = sdist, 
                    hdist = hdist, 
                    link = link, 
                    link.fun = link.fun, 
                    link.inv = link.inv, 
                    link.inv.der = link.inv.der,
                    start = start,
                    start.eta = start.eta,
                    use.optim = use.optim,
                    closed.mle = closed.mle
                    )
}




###### dist_list for normal distribution
if(FALSE) {
  
  dist_list_normal <- list()
  
  par.names <- c("mu", "sigma")
  eta.names <- c("mu", "log(sigma)")
  
  use.optim <- FALSE
  
  
  closed.mle <- function(y, weights = NULL){
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, ny)
    y.w <- (y*weights)[weights != 0]     # just for the calculations where the index of the observations is not of any importance, but zero would have an impact on the result
    ny.w <- length(y.w)
    mu <- mean(y.w)
    sigma <- sqrt(1/ny.w * sum((y.w - mu)^2))
    par <- c(mu, sigma)
    names(par) <- par.names
    return(par)
  }
  
  
  ddist <-  function(y, par, log = TRUE, type = "parameter") {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    if(type == "link") par[2] <- exp(par[2])
    val <- 1/sqrt(2*pi*par[2]^2) * exp(- (y-par[1])^2 / (2*par[2]^2))
    if(log) val <- log(val)
    # val <- dnorm(y, mean = par[1], sd = par[2], log = log)
    return(val)
  }
  
  
  sdist <- function(y, par, type = "parameter") {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    if(type == "parameter") {
      score <- cbind(1/par[2]^2 * (y-par[1]), (-1/par[2] + ((y - par[1])^2)/(par[2]^3)))
      score <- as.matrix(score)
      colnames(score) <- par.names
    }
    if(type == "link") {
      eta <- par
      par <- c(eta[1], exp(eta[2]))                           
      score <- cbind(1/par[2]^2 * (y-par[1]), (-1/par[2] + ((y - par[1])^2)/(par[2]^3)) * exp(eta[2]))
      score <- as.matrix(score)
      colnames(score) <- eta.names
    }
    return(score)
  }
  
  
  hdist <- function(y, par, type = "parameter", weights = NULL) {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, length(y))
    
    if(type == "parameter") {
      d2ldm2 <- sum(weights * rep.int(-1/par[2]^2, ny))
      d2ldmdd <- sum(weights * (-2)*(y-par[1])/par[2]^3)                      # should be 0 for exact parameters (here ~ e-17 due to calculations)
      d2ldd2 <- sum(weights * (1/par[2]^2 - 3*(y-par[1])^2/par[2]^4))         # should be -2/sigma^2 for exact parameters 
      
      hess <- matrix(c(d2ldm2, d2ldmdd, d2ldmdd, d2ldd2), nrow = 2)
      colnames(hess) <- rownames(hess) <-  par.names
    }
    
    if(type == "link") {
      eta <- par
      par <- c(eta[1], exp(eta[2]))                           
      
      d2ld.etamu2 <- sum(weights * rep.int(-1/par[2]^2, ny))
      d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-par[1])/par[2]^2)          # should be 0 for exact parameters (here ~ e-17 due to calculations)
      d2ld.etasigma2 <- sum(weights * (-2)*(y-par[1])^2/par[2]^2)         
      
      hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
      colnames(hess) <- rownames(hess) <-  eta.names
    }
    
    return(hess)
  }
  
  link <- c("identity", "log")
  
  link.fun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- eta.names
    return(eta)
  }
  
  
  link.inv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- par.names
    return(par)
  }
  
  
  link.inv.der <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- names.par
    return(dpardeta)
  }
  
  
  start <- function(y) closed.mle(y)
  
  start.eta <- function(y) {
    start <- closed.mle(y)
    starteta <- c(start[1], log(start[2])) 
    names(starteta) <- eta.names
    return(starteta)
  }
  
  
  dist_list_normal <- list(family.name = "Normal Distribution",
                           ddist = ddist, 
                           sdist = sdist, 
                           hdist = hdist, 
                           link = link, 
                           link.fun = link.fun, 
                           link.inv = link.inv, 
                           link.inv.der = link.inv.der,
                           start = start,
                           start.eta = start.eta,
                           use.optim = use.optim,
                           closed.mle = closed.mle
  )
}







###### dist_list for Poisson distribution
if(FALSE) {
  
  dist_list_poisson <- list()
  
  par.names <- c("mu")
  eta.names <- c("log(mu)")
  
  use.optim <- FALSE
  
  closed.mle <- function(y, weights = NULL){
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, ny)
    y.w <- (y*weights)[weights != 0]     # just for the calculations where the index of the observations is not of any importance, but zero would have an impact on the result
    ny.w <- length(y.w)
    mu <- mean(y.w)
    par <- c(mu)
    names(par) <- par.names
    return(par)
  }
  
  
  ddist <-  function(y, par, log = TRUE, type = "parameter") {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    if(type == "link") par <- exp(par)
    val <- par^y / gamma(y+1) * exp(-par)
    if(log) val <- log(val)
    return(val)
  }
  
  
  sdist <- function(y, par, type = "parameter") {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    if(type == "parameter") {
      score <- y/par - 1
      score <- as.matrix(score)
      colnames(score) <- par.names
    }
    if(type == "link") {
      eta <- par
      par <- exp(eta)                           
      score <- (y - par)
      score <- as.matrix(score)
      colnames(score) <- eta.names
    }
    return(score)
  }
  
  
  hdist <- function(y, par, type = "parameter", weights = NULL) {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, length(y))
    
    if(type == "parameter") {
      hess <- -y/(par^2)
      hess <- as.matrix(sum(weights * hess))
      colnames(hess) <- rownames(hess) <-  par.names
    }
    
    if(type == "link") {
      eta <- par
      par <- exp(eta)                           
      hess <- rep(-par, ny)
      hess <- as.matrix(sum(weights * hess))
      colnames(hess) <- rownames(hess) <-  eta.names
    }
    
    return(hess)
  }
  
  link <- c("log")
  
  link.fun <- function(par) {
    eta <- log(par)
    names(eta) <- eta.names
    return(eta)
  }
  
  
  link.inv <- function(eta) {
    par <- exp(eta)
    names(par) <- par.names
    return(par)
  }
  
  
  link.inv.der <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- names.par
    return(dpardeta)
  }
  
  
  start <- function(y) closed.mle(y)
  
  start.eta <- function(y) {
    start <- closed.mle(y)
    starteta <- log(start) 
    names(starteta) <- eta.names
    return(starteta)
  }
  
  use.optim <- FALSE
  
  
  dist_list_poisson <- list(family.name = "Poisson Distribution",
                           ddist = ddist, 
                           sdist = sdist, 
                           hdist = hdist, 
                           link = link, 
                           link.fun = link.fun, 
                           link.inv = link.inv, 
                           link.inv.der = link.inv.der,
                           start = start,
                           start.eta = start.eta,
                           use.optim = use.optim,
                           closed.mle = closed.mle
  )
}






###### dist_list for Poisson distribution
if(FALSE) {
  
  dist_list_exp <- list()
  
  par.names <- c("lambda")
  eta.names <- c("log(lambda)")
  
  use.optim <- FALSE
  
  closed.mle <- function(y, weights = NULL){
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, ny)
    y.w <- (y*weights)[weights != 0]     # just for the calculations where the index of the observations is not of any importance, but zero would have an impact on the result
    ny.w <- length(y.w)
    lambda <- ny.w / sum(y.w)
    par <- c(lambda)
    names(par) <- par.names
    return(par)
  }
  
  
  ddist <-  function(y, par, log = TRUE, type = "parameter") {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    if(type == "link") par <- exp(par)
    val <- par * exp(-par * y)
    if(log) val <- log(val)
    return(val)
  }
  
  
  sdist <- function(y, par, type = "parameter") {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    if(type == "parameter") {
      score <- 1/par - y
      score <- as.matrix(score)
      colnames(score) <- par.names
    }
    if(type == "link") {
      eta <- par
      par <- exp(eta)                           
      score <- 1 - y * par
      score <- as.matrix(score)
      colnames(score) <- eta.names
    }
    return(score)
  }
  
  
  hdist <- function(y, par, type = "parameter", weights = NULL) {     # input parameter have to fit the type (par or eta) but in the input line denoted by par for both cases
    ny <- length(y)
    if(is.null(weights)) weights <- rep.int(1, length(y))
    
    if(type == "parameter") {
      hess <- -1/(par^2)
      hess <- as.matrix(sum(weights * hess))
      colnames(hess) <- rownames(hess) <-  par.names
    }
    
    if(type == "link") {
      eta <- par
      par <- exp(eta)                           
      hess <- -y * par
      hess <- as.matrix(sum(weights * hess))
      colnames(hess) <- rownames(hess) <-  eta.names
    }
    
    return(hess)
  }
  
  link <- c("log")
  
  link.fun <- function(par) {
    eta <- log(par)
    names(eta) <- eta.names
    return(eta)
  }
  
  
  link.inv <- function(eta) {
    par <- exp(eta)
    names(par) <- par.names
    return(par)
  }
  
  
  link.inv.der <- function(eta) {
    dpardeta <- exp(eta)
    names(dpardeta) <- names.par
    return(dpardeta)
  }
  
  
  start <- function(y) closed.mle(y)
  
  start.eta <- function(y) {
    start <- closed.mle(y)
    starteta <- log(start) 
    names(starteta) <- eta.names
    return(starteta)
  }
  
  use.optim <- FALSE
  
  
  dist_list_exp <- list(family.name = "Exponential Distribution",
                            ddist = ddist, 
                            sdist = sdist, 
                            hdist = hdist, 
                            link = link, 
                            link.fun = link.fun, 
                            link.inv = link.inv, 
                            link.inv.der = link.inv.der,
                            start = start,
                            start.eta = start.eta,
                            use.optim = use.optim,
                            closed.mle = closed.mle
  )
}