absc <- c(1e-10,.1,.15,.21,.25,.4,.5,1 - 1e-10)

n <- 100; x <- 20;

fbin <- function(p){
  dbinom(x,n,p,log = T)
}

f <- sapply(absc,fbin)

fmxind <- which(f==max(f))

if (f[fmxind] > fbin(absc[fmxind] - 1e-5) & f[fmxind] > fbin(absc[fmxind] + 1e-5)) {
  maxind <- 0
} else if (f[fmxind] < fbin(absc[fmxind] + 1e-5)) {
  maxind <- 1
} else if (fbin(absc[fmxind] - 1e-5) > f[fmxind]) {
  maxind <- 2
}

### calculating intersection x for middle section (containing max); using Rcpp function intsct2
if (maxind==2) {
  zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
} else {
  zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
}

zval <- intsct2(zabsc,x,n)

abscaug <- sort(unique(c(absc,zval)))
# f3mxindaug <- unique(c(f3mxind,f3mxind + 1))
### "chord" function for lower hull
flowerhull <- function(p,absc,f) {
  selind <- tail(which(absc <= p),1)
  xcap <- p - absc[selind]
  exp(f[selind] + (f[selind + 1] - f[selind])/(absc[selind + 1] - absc[selind]) * xcap)
}

f <- sapply(abscaug,fbin)
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
###               Upper hull
#########################################################################################
#########################################################################################
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
  } else if ((abscaug[k]==zval & zind==1) | (abscaug[k + 1]==zval & zind==1 & maxind==0) | 
             (abscaug[k]==zval & zind==0 & maxind==0) | (abscaug[k] > zval)) {
    fval <- exp(f[k])
  } 
  fval
}

#########################################################################################
#########################################################################################
###     Setting up the probabilities per "segment"
#########################################################################################
#########################################################################################
fliksum <- function(abscaug,f,zval) {
  lvec <- NULL
  
  zind <- ifelse(fbin(zval + 1e-5) > fbin(zval), 0,1) ### indicator 1 if maximum after zval
  
  for (k in seq_along(abscaug[-1])) {
    
    if (k%in%f3mxind) {
      if ((abscaug[k + 1] < zval) | (abscaug[k + 1] == zval & zind==0)) {
        prob <- exp(f[k + 1])*(abscaug[k + 1] - abscaug[k])
      } else if (abscaug[k] < zval & zind==1) {
        f0 <- f[k - 1]
        f1 <- f[k]
        p0 <- abscaug[k - 1]
        p1 <- abscaug[k]
        a <- (f1 - f0)/(p1 - p0)
        fval <- exp(f1)*(exp(a*(zval - p1)) - 1)/a
      } else if ((abscaug[k] == zval & zind==0) | (abscaug[k] > zval)) {
        prob <- exp(f[k])*(abscaug[k + 1] - abscaug[k])
      } else if (abscaug[k] == zval & zind==0) {
        f0 <- f[k + 1]
        f1 <- f[k + 2]
        p0 <- abscaug[k + 1]
        p1 <- abscaug[k + 2]
        a <- (f1 - f0)/(p1 - p0)
        fval <- exp(f1)*(exp(a*(zval - p1)) - 1)
      }
    } else if (k%in%lowerind) {
      prob <- exp(f[k + 1])*(abscaug[k + 1] - abscaug[k])
    } else if (k%in%upperind) {
      prob <- exp(f[k])*(abscaug[k + 1] - abscaug[k])
    } 
    lvec[k] <- prob
  }
  lvec
}
### CDF
lvec <- fliksum(abscaug,f,zval)
#########################################################################################
###    Sampling a p from cummulative distribution
#########################################################################################
fpsample <- function(abscaug,lvec,f,zval) {
  zind <- ifelse(fbin(zval + 1e-5) > fbin(zval), 0,1) ### indicator 1 if maximum after zval

  pvec <- lvec/sum(lvec)
  pveccum <- c(0,cumsum(pvec))
  
  pran <- runif(1)
  k <- max(which(pveccum <= pran))
  pbase <- pveccum[k - 1]
  ptop <- pveccum[k]
  prm <- pran - pbase
  
  if ((abscaug[k + 1] == zval & zind==0) | (abscaug[k] > zval) | (abscaug[k] == zval & zind==1) | (abscaug[k + 1] < zval)) {
    pout <- abscaug[k] + (abscaug[k + 1] - abscaug[k])*prm/(ptop - pbase)
  } else if (abscaug[k + 1] == zval & zind==1) {
    f0 <- f[k - 1]
    f1 <- f[k]
    p0 <- abscaug[k - 1]
    p1 <- abscaug[k]
    a <- (f1 - f0)/(p1 - p0)
    pout <- (log(exp(f1) + a*prm*sum(lvec)) + a*p1 - f1)/a
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
