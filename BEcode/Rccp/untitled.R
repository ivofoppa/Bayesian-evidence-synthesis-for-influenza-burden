fprob <- function(abscaug,f,zval) {
  pvec <- NULL
  
  zind <- ifelse(fbin(zval + 1e-5) > fbin(zval), 0,1) ### indicator 1 if maximum after zval
  
  for (k in seq_along(abscaug[-1])) {
 
        if (abscaug[k] < zval & zind==0) {
        prob <- exp(f[k + 1])*(abscaug[k + 1] - abscaug[k])
        } else if (abscaug[k] > zval){
          prob <- exp(f[k])*(abscaug[k + 1] - abscaug[k])
        } else if (abscaug[k + 1] == zval & zind==1){
          f0 <- f[k - 1]
          f1 <- f[k]
          p0 <- abscaug[k - 1]
          p1 <- abscaug[k]
          a <- (f1 - f0)/(p1 - p0)
          prob <- exp(f1)*(exp(a*(zval - p1)) - 1)/a
        } else if (abscaug[k] == zval & zind==0){
          f0 <- f[k + 1]
          f1 <- f[k + 2]
          p0 <- abscaug[k + 1]
          p1 <- abscaug[k + 2]
          a <- (f1 - f0)/(p1 - p0)
          prob <- exp(f0)*(1 - exp(a*(zval - p0)))/a
        } else if (abscaug[k] > zval & zind==1){
          prob <- exp(f[k])*(abscaug[k + 1] - abscaug[k])
        }
    pvec[k] <- prob
  }
  pvec
}

fpsample <- function(abscaug,lvec,f,zval) {
  pvec <- lvec/sum(lvec)
  pveccum <- cumsum(pvec)
  pran <- runif(1)
  k <- max(which(c(0,pveccum) <= pran))
  pbase <- pveccum[k - 1]
  ptop <- pveccum[k]
  prm <- pran - pbase
  
  zind <- ifelse(fbin(zval + 1e-5) > fbin(zval), 0,1) ### indicator 1 if maximum after zval
  
  if (abscaug[k + 1] < zval & zind==0 | abscaug[k] > zval | abscaug[k] == zval & zind==1 ) {
    pout <- abscaug[k] + (abscaug[k + 1] - abscaug[k])*prm/(ptop - pbase)
  } else if (abscaug[k + 1] == zval & zind==1) {
    f0 <- f[k - 1]
    f1 <- f[k]
    p0 <- abscaug[k - 1]
    p1 <- abscaug[k]
    a <- (f1 - f0)/(p1 - p0)
    pout <- (log(a + a*prm*sum(lvec)) + a*p1 - f1)/a
  } else if (abscaug[k] == zval & zind==0) {
    f0 <- f[k + 1]
    f1 <- f[k + 2]
    p0 <- abscaug[k + 1]
    p1 <- abscaug[k + 2]
    a <- (f1 - f0)/(p1 - p0)
    pout <- (log(exp(f0 + a*(zval - p0)) + prm*a*sum(lvec)) + a*p0 - f0)/a
  }
  pout
}

