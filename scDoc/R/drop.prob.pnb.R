#' Calculate dropout probability for PNB model
#' 
#' @description Function to calculate dropout probability 
#' 
#' @param pars vector containing parameters, gamma, alpha, and delta
#' @param s PNB mixture parameter
#' @param y vector containing observed expression values
#' @param offset vector containing log effective library size
#' @param max.val maximum value for computation
#' @param tol tolerance for computation
#' 
#' @return a vector containing dropout probability
#' 
#' @export
#' 



drop.prob.pnb <- function(pars, s, y, offset, max.val=1e100, tol=1e-10) {
  
  gamma <- pars[1]
  alpha <- pars[2]
  delta <- pars[3]
  
  lambda <- pmin(exp(gamma + offset), max.val)
  mu <- pmin(exp(alpha + offset), max.val)
  phi <- max(min(exp(delta), max.val), tol)
  
  d <- (s*dpois(y,lambda))/(s*dpois(y,lambda) + (1-s)*dnbinom(y,mu=mu,size=phi))
  
  return(d)
}