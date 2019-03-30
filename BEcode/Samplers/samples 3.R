setwd('C:/Users/vor1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/BEcode/Samplers')
source('binom derivative free.R')
# source('binom derivative free beta.R')
#########################################################################################
#########################################################################################
###  Drawing samples from posterior    ##################################################
#########################################################################################
#########################################################################################
n <- 100
x <- 1
nsim <- 10000
p0 <- .2

absc <- c(0,0.02,.05,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9,1)

absc <- abscPrep(absc,x,n)
zabsc <- zabsc(absc,x,n)
zval <- intsct2(zabsc,x,n)

fbin <- function(p){
  dbinom(x,n,p,log = T)
}

absc <- unique(sort(c(absc,zval)))

f <- sapply(absc,fbin)

fmxind <- which(f==max(f))

lvec <- fliksum(absc,f,zval,maxind)

probls <- NULL 
crit <- .75 

while (length(probls) < nsim) {
  psample <- fpsample(abscaug,lvec,f,zval,maxind)
  aran <- runif(1)
  lhull <- flowerhull(psample,abscaug,f)
  uhull <- fupperhull(psample,abscaug,f,zval,maxind)
  if (aran < lhull/uhull) {
    probls <- c(probls,psample)
  } else {
    aratio <- exp(fbin(psample))/uhull
    if (aran < aratio) {
      probls <- c(probls,psample)
    } else if (aratio < crit ) {
      abscaug <- unique(sort(round(c(abscaug,psample),digits = 8)))
      f <- sapply(abscaug,fbin)
      
      fmxind <- which(f==max(f))
      
      if (f[fmxind] > fbin(max(abscaug[fmxind] - 1e-5,0)) & f[fmxind] > fbin(min(abscaug[fmxind] + 1e-5,1))) {
        maxind <- 0
      } else if (f[fmxind] < fbin(min(abscaug[fmxind] + 1e-5,1))) {
        maxind <- 1
      } else if (fbin(max(abscaug[fmxind] - 1e-5,0)) > f[fmxind]) {
        maxind <- 2
      }
      
      ### calculating intersection x for middle section (containing max); using Rcpp function intsct2
      if (maxind==2) {
          zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
        } else {
          zabsc <- abscaug[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
        } 
      
      zval <- intsct2(zabsc,x,n)
      
      abscaug <- sort(unique(c(abscaug,zval)))
      f <- sapply(abscaug,fbin)
      lvec <- fliksum(abscaug,f,zval,maxind)
    }
  }
}

mean(probls)
