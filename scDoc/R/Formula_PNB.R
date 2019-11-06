#=======================
# All Formula for PNB
#=======================

# -1 * loglik for optim function

# loglikelihood function - Poisson-NB distribution
loglik.PNB <- function(pars, s, y, offset, max.val=1e100, tol=1e-10) {
  
  gamma <- pars[1]
  alpha <- pars[2]
  delta <- pars[3]
  
  lambda <- pmin(exp(gamma + offset), max.val)
  mu <- pmin(exp(alpha + offset), max.val)
  phi <- max(min(exp(delta), max.val), tol)
  
  loglik <- sum(log(s*dpois(y,lambda) + (1-s)*dnbinom(y,mu=mu,size=phi)))
  
  return(-1*loglik)
}


# loglikelihood function - joint model
loglik.PNB.complete <- function(pars, s, z, y, offset, max.val=1e100, tol=1e-10) {
  
  gamma <- pars[1]
  alpha <- pars[2]
  delta <- pars[3]
  
  lambda <- pmin(exp(gamma + offset), max.val)
  mu <- pmin(exp(alpha + offset), max.val)
  phi <- max(min(exp(delta), max.val), tol)
  
  l1c <- sum(z*log(s)+(1-z)*log(1-s))
  l2c <- sum(z*dpois(y, lambda, log=T) + (1-z)*dnbinom(y, mu=mu, size=phi, log=T))
  
  return(-1*(l1c+l2c))
}


# gradient - joint model
loglik.PNB.grad <- function(pars, s, z, y, offset, max.val=1e100, tol=1e-10) {
  
  gamma <- pars[1]
  alpha <- pars[2]
  delta <- pars[3]
  
  lambda <- pmin(exp(gamma + offset), max.val)
  mu <- pmin(exp(alpha + offset), max.val)
  phi <- max(min(exp(delta), max.val), tol)
  
  grad.gamma <- sum(z*(y-lambda))
  grad.alpha <- sum((1-z)*(y-(y+phi)*mu/(mu+phi)))
  grad.delta <- sum((1-z)*(digamma(phi+y) - digamma(phi) + 1 + log(phi) - (y+phi)/(mu+phi) - log(mu+phi)))*phi
  
  return(-1*c(grad.gamma,
              grad.alpha,
              grad.delta))
}

