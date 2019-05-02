mu <- 0.6
sd <- .1

fit_beta <- function(mu,sd){
  beta0 <- 1
  alpha0 <- mu*beta0/(1 - mu)
  var0 <- sqrt(alpha0*beta0/((alpha0 + beta0)^2 * (alpha0 + beta0 + 1)))
}

