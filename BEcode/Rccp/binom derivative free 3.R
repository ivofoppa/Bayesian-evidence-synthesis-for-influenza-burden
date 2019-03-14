absc <- c(1e-10,.1,.15,.2,.25,.4,.5,1 - 1e-10)

n <- 100; x <- 20;

fbin <- function(p){
  dbinom(x,n,p,log = T)
}

f <- sapply(absc,fbin)

f3mx <- sort(f)[(length(f) - 2):length(f)]

f3mxind <- which(f%in%f3mx)

if (fbin(absc[f3mxind[2]] + 1e-5) > fbin(absc[f3mxind[2]])) {
  f3mxind <- f3mxind[-1]
} else f3mxind <- f3mxind[-3]

f3mxindplus <- unique(c(f3mxind - 1,f3mxind + 1))
absc2 <- absc[f3mxindplus]

### calculating intersection x for middle section (containing max); using Rcpp function intsct2
zval <- intsct2(absc2,x,n)

### Indices of lower and upper ranges where "step method" can be used
lowerind <- seq(1,min(f3mxind) - 1)
### Augmenting absc for use in sampling
abscaug <- unique(c(absc[lowerind],absc[f3mxind[1]],zval,absc[f3mxind[2]],absc[upperind - 1],max(absc)))
f3mxindaug <- unique(c(f3mxind,f3mxind + 1))
upperind <- seq(max(f3mxindaug) + 1,length(abscaug))

f <- sapply(abscaug,fbin)

### Determining location of max and calculating middle prob, using "chord method"
if (fbin(zval) > fbin(zval + 1e-5)) {
  k <- f3mxindaug[1]
  a <- (f[k] - f[k - 1])/(absc[k] - absc[k - 1]) 
  pvecmiddle <- exp(f[k])/a * (exp(a * (zval - absc[k])) - 1)
  ### Indices for middle probs
} else if (fbin(zval) < fbin(zval + 1e-5)) {
  k <- f3mxind[2]
  b <- (f[k + 1] - f[k - 1])/(abscaug[k + 1] - abscaug[k - 1]) 
  pvecmiddle <- exp(f[k])/b * (1 - exp(b * (zval - absc[k + 1])))
  ### Indices for middle probs
}
### Udating f
### Probs of lower indices, using "step method"
pveclower <- NULL
for (k in lowerind){
  plower <- exp(f[k + 1]) * (abscaug[k + 1] - abscaug[k])
  pveclower <- c(pveclower,plower)
}
### Probs of upper indices, using "step method"
pvecupper <- NULL
for (k in head(upperind,-1)){
  pupper <- exp(f[k]) * (abscaug[k] - abscaug[k + 1])
  pvecupper <- c(pvecupper,pupper)
}

if (fbin(zval + 1e-5) > fbin(zval)) {
  k <- f3mxind[1]
  plower <- exp(f[k]) * (abscaug[k] - abscaug[k - 1])
  pveclower <- c(pveclower,plower)
} else {
  k <- f3mxind[2]
  plower <- exp(f[k]) * (abscaug[k] - abscaug[k - 1])
  pveclower <- c(pveclower)
}

### putting the values together and normalizing
pvec <- c(pveclower, pvecmiddle, pvecupper)
pnorm <- sum(pvec)
pvec <- pvec /pnorm
pveccum <- cumsum(pvec)

### "chord" function for lower hull
flowerhull <- function(p,absc,f) {
  selind <- tail(which(absc <= p),1)
  xcap <- p - absc[selind]
  exp(f[selind] + (f[selind + 1] - f[selind])/(absc[selind + 1] - absc[selind]) * xcap)
}

f <- sapply(abscaug,fbin)
#########################################################################################
#########################################################################################
###               Upper hull
#########################################################################################
#########################################################################################
fupperhull <- function(p,abscaug,f) {
  selind <- max(which(abscaug <= p))
  
if (selind%in%f3mxind) {
    if (p <= zval) {
      f0 <- f[selind - 1]
      f1 <- f[selind]
      p0 <- abscaug[selind - 1]
      p1 <- abscaug[selind]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f1 + a*(p - p1))
    } else if (p > zval) {
      f0 <- f[selind + 1]
      f1 <- f[selind + 2]
      p0 <- abscaug[selind + 1]
      p1 <- abscaug[selind + 2]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f0 + a*(p - p0))
    }
  } else if (selind%in%lowerind) {
    fexp <- exp(f[selind + 1])
    fval <- fexp
    
  } else if (selind%in%upperind) {
    fexp <- exp(f[selind ])
    fval <- fexp
  }
  fval
}

### Sampling
rval <- runif(1)
pselind <- max(which(pveccum <= rval))
abscsel <- abscaug[pselind]

prm <- rval - pveccum[pselind]

if (pselind==f3mxind[1]) {
  f0 <- f[pselind - 1]
  f1 <- f[pselind]
  p0 <- abscaug[pselind - 1]
  p1 <- abscaug[pselind]
  a <- (f1 - f0)/(p1 - p0)
  p <- (log(prm*pnorm*a + exp(f1)) + a*p1 - f1)/a
} else if (pselind==f3mxind[2]) {
  f0 <- f[pselind]
  f1 <- f[pselind + 1]
  p0 <- abscaug[pselind]
  p1 <- abscaug[pselind + 1]
  a <- (f1 - f0)/(p1 - p0)
  p <- (log(exp(f0 + a*(zval - po)) - prm*pnorm*a) + a*x0 - f1)/a
} else {
  
}

