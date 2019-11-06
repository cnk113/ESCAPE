#' Possion-NB Mixture Model Parameter Estimate
#'
#' @description Function to estimate parameters for Poisson-NB
#'              mixutre model. It uses "BFGS", a quasi-Newton 
#'              method and depends on optim function.
#'              
#' @param y vector containing expression value for a gene
#' @param offsets vector containing log transformed effective 
#'                library size. It could be total counts per 
#'                sample or calculated effective library.
#' @param ini.pars list containing initial parameters 
#'                 provided by user. If any of them is undefined, 
#'                 program will estimate initial parameters from 
#'                 data.
#' @param max.val maximum value for computation
#' @param tol tolerance (minimum value) for computation
#' @param EM.maxIter maximum iteration for EM algorithm
#' @param EM.tol convergence threshold for EM algorithm
#' @param ... other arguments passed to optim function. See help 
#'            for optim function.   
#'            
#' @return a list containing estimate results:
#' \itemize{
#'   \item gamma estimation for gamma
#'   \item alpha estimation for alpha
#'   \item delta estimation for delta
#'   \item s estimation for s
#'   \item z estimation for z
#'   \item loglik log-likehood
#'   \item EM.convergence 0, EM algorithm is converged;
#'                        1, not converged.
#'   \item iter number of iteration used in EM algorithm
#' }  
#' 
#' @examples 
#' 
#' # generate dataset from Piosson-NB mixture model
#' pois.nb.mix <- function(n, pars, prob=0.2, offsets, max_val=1e100, tol=1e-10) {
#' 
#'   gamma <- pars[1]
#'   alpha <- pars[2]
#'   delta <- pars[3]
#' 
#'   lambda <- pmin(exp(gamma + offsets - mean(offsets)), max_val)
#'   mu <- pmin(exp(alpha + offsets- mean(offsets)), max_val)
#'   phi <- max(min(exp(delta), max_val), tol)
#' 
#'   u <- runif(n) 
#'   out <- apply( as.matrix(u), 1, function(x) ifelse(x<=prob, rpois(1, lambda), rnbinom(1, mu=mu, size=phi) ) )
#'   out 
#' }
#' 
#' # generate a test sample
#' set.seed(2018)
#' test <- pois.nb.mix(1000, c(-0.5, 7, 1), 0.2, offsets = rep(0,1000))
#' 
#' rs <- param.est(y=test, offsets = rep(0, 1000))
#' rs$gamma
#' rs$alpha
#' rs$delta
#' rs$s
#' rs$iter
#' rs$EM.convergence
#' 
#' @export
#' 
#'               
#'                                           

param.est.pnb <- function(y, offsets, ini.pars=list(gamma=NULL, alpha=NULL, delta=NULL, s=NULL), 
                      max.val=1e100, tol=1e-10, EM.maxIter=300, EM.tol=1e-04, ...) {
  
  require(MASS)
  
  # initial values
  if (!any(unlist(lapply(ini.pars, is.null)))) {
    
    pars <- c(ini.pars$gamma, ini.pars$alpha, ini.pars$delta)
    s.t <- ini.pars$s
    
  } else {
  
    v <- 2
    
    outb <- y <= v
    
    if (all(offsets == 0)) {
      
      if (all(outb == F)) {
        
        ini.s <- 0.5
        
      } else {
        
        fit1 <- suppressWarnings(glm(outb ~ 1, family = "binomial"))
        beta0 <- as.numeric(coef(fit1))
        ini.s <- 1/(1 + exp(-1*beta0))
        
      }
      
    } else {
      
      if (all(outb == F)) {
        
        beta0 <- beta1 <- 0
        
      } else {
        
        fit2 <- suppressWarnings(glm(outb ~ offsets, family = "binomial"))
        beta <- as.numeric(coef(fit2))
        beta0 <- beta[1]
        beta1 <- beta[2]
        
        if (is.na(beta1)) {
          beta1 <- 0
        }
        
      }
      
      ini.s <- 1/(1 + exp(-1*beta0 - beta1*offsets))
      
    }
    
    yp <- y[y<=v]
    ynb <- y[y>v]
    
    if (length(yp) > 0) {
      
      ini.pois <- max(suppressWarnings(fitdistr(yp, "poisson")$estimate), 1e-4)
      
    } else {
      
      ini.pois <- 1e-2
      
    }
    
    znb.parms <- suppressWarnings(try(fitdistr(ynb,"negative binomial"), silent = T))
    
    if(!inherits(znb.parms, "try-error")) {
      
      parms <- znb.parms$estimate
      
      ini.mu <- parms[2]
      ini.phi <- parms[1]
      
    } else {
      ini.mu <- mean(y[y>0])
      ini.phi <- max((var(y[y>0]) - mean(y[y>0]))/ (mean(y[y>0]))^2, tol)
    }
    
    ini.gamma <- mean(log(ini.pois) - offsets)
    ini.alpha <- mean(log(ini.mu) - offsets)
    ini.delta <- max(as.numeric(log(ini.phi)), 1+tol)
    
    pars <- c(ini.gamma, ini.alpha, ini.delta)
    s.t <- ini.s
  }
  
  # E-M
  status <- 0
  loglik <- -Inf
  z.t <- iter <- EM.convergence <- NA
  
  for (iter in 1:EM.maxIter) {
    loglik_pre <- loglik
    
    # E-step
    lambda.t <- pmin(exp(pars[1] + offsets), max.val)
    mu.t <- pmin(exp(pars[2] + offsets), max.val)
    phi.t <- max(min(exp(pars[3]), max.val), tol)
    z.t <- s.t * dpois(y, lambda.t) / (s.t * dpois(y, lambda.t) + (1-s.t) * dnbinom(y, mu=mu.t, size=phi.t))
    z.t[is.na(z.t)] <- 0.5
    
    # M-step
    out.t <- z.t > 0.5
    
    if (all(offsets==0)) {
      
      if (all(out.t == F)) {
        
        s.t <- s.t
        
      } else {
        
        fit3 <- suppressWarnings(glm(out.t ~ 1, family = "binomial"))
        beta0 <- as.numeric(coef(fit3))
        s.t <- 1/(1 + exp(-1*beta0))
        
      }
      
    } else {
      
      if (all(out.t == F)) {
        
        s.t <- s.t
        
      } else {
        
        fit4 <- suppressWarnings(glm(out.t ~ offsets, family = "binomial"))
        beta <- as.numeric(coef(fit4))
        beta0 <- beta[1]
        beta1 <- beta[2]
        
        if (is.na(beta1)) {
          beta1 <- 0
        }
        
        s.t <- 1/(1 + exp(-1*beta0 - beta1*offsets))
        
      }
    }
    
    obj <- try(optim(par=pars,
                     f=loglik.PNB.complete,
                     gr=loglik.PNB.grad,
                     y=y,
                     offset=offsets,
                     s=s.t,
                     z=z.t,
                     method = "BFGS"), 
               silent = T)
    
    
    if(!inherits(obj,'try-error')) {
      
      pars <- obj$par
      loglik <- -1*loglik.PNB(pars=pars, s=s.t, y=y, offset=offsets)
      diff <- (abs(loglik-loglik_pre) < EM.tol)
      
    } else {
      
      diff <- FALSE
      
    }
    
    if (diff) {
      status <- 1
      break
    }
    
  }
  
  EM.convergence <- 1-status
  
  res <- list(gamma=pars[1],
              alpha=pars[2],
              delta=pars[3],
              s=s.t,
              z=z.t,
              loglik=loglik,
              EM.convergence=EM.convergence,
              iter=iter)
  
  return(res)
}