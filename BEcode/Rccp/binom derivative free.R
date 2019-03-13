absc <- c(1e-10,.1,.15,.2,.25,.4,.5,1 - 1e-10)

n <- 100; x <- 20;

fbin <- function(p){
  dbinom(x,n,p,log = T)
}

f <- sapply(absc,fbin)

f3mx <- sort(f)[(length(f) - 2):length(f)]

f3mxind <- which(f%in%f3mx)

if (dbinom(x,n,absc[f3mxind[2]] - 1e-10,log = T) < f[f3mxind[2]]) {
  f3mxind <- f3mxind[-1]
} else f3mxind <- f3mxind[-3]

f3mxindplus <- unique(c(f3mxind - 1,f3mxind + 1))
f3mxindmin <- head(f3mxind,-1) 

absc2 <- absc[f3mxindplus]

### calculating intersection x for middle section (containing max); using Rcpp function intsct2
zval <- intsct2(absc2,x,n)
### calculating probs x_k-int, int-x_(k+1), using "elongated chords" of log lik
pvecmiddle <- NULL
k <- f3mxind[1]
a <- (f[k] - f[k - 1])/(absc[k] - absc[k - 1]) 

b <- (f[k + 2] - f[k + 1])/(absc[k + 2] - absc[k + 1]) 

pvecmiddle[1] <- exp(f[k])/a * (exp(a * (zval - absc[k])) - 1)
pvecmiddle[2] <- exp(f[k + 1])/b * (1 - exp(b * (zval - absc[k + 1])))

### Indices of lower and upper ranges where "step method" can be used
lowerind <- seq(1,min(f3mxind) - 1)
upperind <- seq(max(f3mxind) + 1,length(absc))
### Augmenting absc for use in sampling
abscaug <- unique(c(absc[lowerind],absc[f3mxind[1]],zval,absc[f3mxind[2]],absc[upperind],max(absc)))
### Probs of lower indices, using "step method"
pveclower <- NULL
for (k in lowerind){
  plower <- exp(f[k + 1]) * (absc[k + 1] - absc[k])
  pveclower <- c(pveclower,plower)
}

### Probs of upper indices, using "step method"
pvecupper <- NULL
for (k in upperind){
  pupper <- exp(f[k - 1]) * (absc[k] - absc[k - 1])
  pvecupper <- c(pvecupper,pupper)
}
### Indices for lower probs
pmiddleind <- f3mxind
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

