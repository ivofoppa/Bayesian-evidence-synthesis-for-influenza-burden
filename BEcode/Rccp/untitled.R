fupperhull <- function(p,abscaug,f,zval,maxind) {
  
  k <- max(which(abscaug < p))
  zind <- ifelse(fbin(zval - 1e-5) > fbin(zval) , 1, 0)
  if ((abscaug[k + 1]==zval & zind==0) | (abscaug[k + 1]<zval )){
    fval <- exp(f[k + 1])
  } else if (abscaug[k + 1] == zval & zind==1 & maxind!=0){
    f0 <- f[k - 1]
    f1 <- f[k]
    p0 <- abscaug[k - 1]
    p1 <- abscaug[k]
    a <- (f1 - f0)/(p1 - p0)
    fval <- exp(f1 + a*(p - p1))
  } else if (abscaug[k] == zval & zind==0 & maxind!=0){
    f0 <- f[k + 1]
    f1 <- f[k + 2]
    p0 <- abscaug[k + 1]
    p1 <- abscaug[k + 2]
    a <- (f1 - f0)/(p1 - p0)
    fval <- exp(f0 + a*(p - p0))
  } else if ((abscaug[k]==zval & zind==1) | (abscaug[k + 1]==zval & zind==1 & maxind==0) | (abscaug[k]==zval & zind==0 & maxind==0) | (abscaug[k] > zval)) {
    fval <- exp(f[k])
  } 
  fval
}

fliksum <- function(abscaug,f,zval) {
  lvec <- NULL
  
  zind <- ifelse(fbin(zval - 1e-5) > fbin(zval) , 1, 0)
  
  for (k in seq_along(abscaug[-1])) {
    
    if ((abscaug[k + 1]==zval & zind==0) | (abscaug[k + 1]<zval )){
      prob <- exp(f[k + 1])*(abscaug[k + 1] - abscaug[k])
    } else if (abscaug[k + 1] == zval & zind==1 & maxind!=0){
      f0 <- f[k - 1]
      f1 <- f[k]
      p0 <- abscaug[k - 1]
      p1 <- abscaug[k]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f1)*(exp(a*(zval - p1)) - 1)/a
    } else if (abscaug[k] == zval & zind==0 & maxind!=0){
      f0 <- f[k + 1]
      f1 <- f[k + 2]
      p0 <- abscaug[k + 1]
      p1 <- abscaug[k + 2]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f1)*(exp(a*(zval - p1)) - 1)
    } else if ((abscaug[k]==zval & zind==1) | (abscaug[k + 1]==zval & zind==1 & maxind==0) | (abscaug[k]==zval & zind==0 & maxind==0) | (abscaug[k] > zval)) {
      prob <- exp(f[k])*(abscaug[k + 1] - abscaug[k])
    } 
    lvec[k] <- prob
  }
  lvec
}
