library(Rcpp)
setwd('C:/Users/vor1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/BEcode/Samplers')

sourceCpp('intsct2.cpp')
x <- 20; p <- .20; 

absc <- c(seq(x,200,10),1000)


fbin <- function(n){
  dbinom(x,n,p,log = T)
}

f <- sapply(absc,fbin)

fmxind <- which(f==max(f))

if (f[fmxind] > fbin(absc[max(fmxind - 1,1)]) & f[fmxind] > fbin(absc[min(fmxind + 1,length(absc))])) {
  maxind <- 0
} else if (f[fmxind] < fbin(absc[min(fmxind + 1,length(absc))])) {
  maxind <- 1
} else if (fbin(absc[max(fmxind - 1,1)]) > f[fmxind]) {
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
#########################################################################################
#########################################################################################
### "chord" function for lower hull
#########################################################################################
#########################################################################################
flowerhulln <- function(n,absc,f) {
  selind <- tail(which(absc <= n),1)
  xcap <- n - absc[selind]
  if (f[selind] > -Inf) {
    exp(f[selind] + (f[selind + 1] - f[selind])/(abscaug[selind + 1] - abscaug[selind]) * xcap)
  } else 0
}

f <- sapply(abscaug,fbin)
#########################################################################################
#########################################################################################
###               Upper hull
#########################################################################################
#########################################################################################
fupperhull <- function(n,abscaug,f,zval,maxind) {
  
  k <- max(which(abscaug < n))
  zind <- ifelse(fbin(max(zval - 1e-5,zval/2)) > fbin(zval) , 1, 0)
  if ((abscaug[k + 1]==zval & zind==0) | (abscaug[k + 1]<zval )){
    fval <- exp(f[k + 1])
  } else if (abscaug[k + 1] == zval & zind==1 & maxind!=0){
    f0 <- f[k - 1]
    f1 <- f[k]
    n0 <- abscaug[k - 1]
    p1 <- abscaug[k]
    a <- (f1 - f0)/(p1 - n0)
    fval <- exp(f1 + a*(n - p1))
  } else if (abscaug[k] == zval & zind==0 & maxind!=0){
    f0 <- f[k + 1]
    f1 <- f[k + 2]
    n0 <- abscaug[k + 1]
    p1 <- abscaug[k + 2]
    a <- (f1 - f0)/(p1 - n0)
    fval <- exp(f0 + a*(n - n0))
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
fliksum <- function(abscaug,f,zval,maxind) {
  lvec <- NULL
  
  zind <- ifelse(fbin(max(zval - 1e-5,zval/2)) > fbin(zval) , 1, 0)
  
  for (k in seq_along(abscaug[-1])) {
    
    if ((abscaug[k + 1]==zval & zind==0) | (abscaug[k + 1]<zval )){
      prob <- exp(f[k + 1])*(abscaug[k + 1] - abscaug[k])
    } else if (abscaug[k + 1] == zval & zind==1 & maxind!=0){
      f0 <- f[k - 1]
      f1 <- f[k]
      p0 <- abscaug[k - 1]
      p1 <- abscaug[k]
      a <- (f1 - f0)/(p1 - p0)
      prob <- exp(f1)*(exp(a*(zval - p1)) - 1)/a
    } else if (abscaug[k] == zval & zind==0 & maxind!=0){
      f0 <- f[k + 1]
      f1 <- f[k + 2]
      p0 <- abscaug[k + 1]
      p1 <- abscaug[k + 2]
      a <- (f1 - f0)/(p1 - p0)
      prob <- exp(f0)*(1 - exp(a*(zval - p0)))/a
    } else if ((abscaug[k]==zval & zind==1) | (abscaug[k + 1]==zval & zind==1 & maxind==0) | 
               (abscaug[k]==zval & zind==0 & maxind==0) | (abscaug[k] > zval)) {
      prob <- exp(f[k])*(abscaug[k + 1] - abscaug[k])
    } 
    lvec[k] <- prob
  }
  lvec
}
### CDF
lvec <- fliksum(abscaug,f,zval,maxind)
#########################################################################################
###    Sampling a p from cummulative distribution
#########################################################################################
fpsample <- function(abscaug,lvec,f,zval,maxind) {
  zind <- ifelse(fbin(max(zval - 1e-5,zval/2)) > fbin(zval) , 1, 0)
  
  pvec <- lvec/sum(lvec)
  pveccum <- c(0,cumsum(pvec))
  
  pran <- runif(1)
  k <- max(which(pveccum <= pran))
  pbase <- pveccum[k]
  ptop <- pveccum[k + 1]
  prm <- pran - pbase
  
  if (abscaug[k + 1] == zval & zind==1 & maxind > 0) {
    f0 <- f[k - 1]
    f1 <- f[k]
    p0 <- abscaug[k - 1]
    p1 <- abscaug[k]
    a <- (f1 - f0)/(p1 - p0)
    pout <- (log(exp(f1) + a*prm*sum(lvec)) + a*p1 - f1)/a
  } else if (abscaug[k] == zval & zind==0 & maxind > 0) {
    f0 <- f[k + 1]
    f1 <- f[k + 2]
    p0 <- abscaug[k + 1]
    p1 <- abscaug[k + 2]
    a <- (f1 - f0)/(p1 - p0)
    pout <- (log(exp(f0 + a*(zval - p0)) + prm*a*sum(lvec)) + a*p0 - f0)/a
  } else {
    pout <- abscaug[k] + (abscaug[k + 1] - abscaug[k])*prm/(ptop - pbase)
  } 
  pout
}
