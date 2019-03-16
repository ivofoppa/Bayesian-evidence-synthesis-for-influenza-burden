fupperhull2 <- function(p,abscaug,f,zval,maxind) {
  k <- max(which(abscaug <= p))
    if ((abscaug[k + 1]==zval & maxind==1) | (abscaug[k + 1]<zval)){
      fval <- exp(f[k + 1])
    } else if (abscaug[k + 1] == zval & maxind==2){
      f0 <- f[k - 1]
      f1 <- f[k]
      p0 <- abscaug[k - 1]
      p1 <- abscaug[k]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f1 + a*(p - p1))
    } else if (abscaug[k] == zval & maxind==1){
      f0 <- f[k + 1]
      f1 <- f[k + 2]
      p0 <- abscaug[k + 1]
      p1 <- abscaug[k + 2]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f0 + a*(p - p0))
    } else if ((abscaug[k]==zval & maxind==2) | (abscaug[k] > zval)){
      fval <- exp(f[k])
    } else if (maxind==0) {
    fval <- ifelse(p <= zval, exp(f[k + 1]),exp(f[k]))
  }
  fval
}
