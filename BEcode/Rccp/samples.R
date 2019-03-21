setwd('C:/Users/vor1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/BEcode/Rccp')
source('binom derivative free 5.R')
# source('binom derivative free beta.R')
#########################################################################################
#########################################################################################
###  Drawing samples from posterior    ##################################################
#########################################################################################
#########################################################################################
n <- 100
xls <- 1:50
nsim <- 100000
p0 <- .2

meanprobls <- NULL

for (x in xls) {
  absc <- c(1e-10,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9,1 - 1e-10)
  

  fbin <- function(p){
    dbinom(x,n,p,log = T)
  }
  f <- sapply(absc,fbin)
  
  fmxind <- which(f==max(f))
  
  if (f[fmxind] > fbin(max(absc[fmxind] - 1e-5,0)) & f[fmxind] > fbin(min(absc[fmxind] + 1e-5,1))) {
    maxind <- 0
  } else if (f[fmxind] < fbin(min(absc[fmxind] + 1e-5,1))) {
    maxind <- 1
  } else if (fbin(max(absc[fmxind] - 1e-5,0)) > f[fmxind]) {
    maxind <- 2
  }
  
  minpexp <- -10
  ### calculating intersection x for middle section (containing max); using Rcpp function intsct2
  if (maxind==2) {
    if ((fmxind > 2 & fmxind < length(absc))){
    zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
    } else if (fmxind <= 2) {
      while (fmxind <= 2) {
        eval(parse(text=paste0('minp <- 1e',minpexp)))
        absc <- sort(unique(c(seq(minp,absc[2],length.out = 4),absc)))
        f <- sapply(absc,fbin)
        fmxind <- which(f==max(f))
        minpexp <- minpexp*2
      }
      zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
    } else if (fmxind==length(absc)) {
      while (fmxind==length(absc)) {
        absc <- sort(unique(c(seq(tail(absc,2)[1],tail(absc,1),length.out = 3),absc)))
      f <- sapply(absc,fbin)
      fmxind <- which(f==max(f))
      }
      zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
    }
  } else {
    if (fmxind > 1 & fmxind < (length(absc) - 1)) {
    zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
    } else if (fmxind == 1) {
      while (fmxind==1) {
        absc <- sort(unique(c(seq(absc[1],absc[2],length.out = 3),absc)))
      f <- sapply(absc,fbin)
      fmxind <- which(f==max(f))
      }
      zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
    } else if (fmxind > (length(absc) - 2)) {
      while (fmxind > (length(absc) - 2)) {
      absc <- sort(unique(c(seq(tail(absc,2)[1],tail(absc,1),length.out = 4),absc)))
      f <- sapply(absc,fbin)
      fmxind <- which(f==max(f))
      }
      zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
    }
  }
  
  zval <- intsct2(zabsc,x,n)
  
  abscaug <- sort(unique(c(absc,zval)))
  
  f <- sapply(abscaug,fbin)
  lvec <- fliksum(abscaug,f,zval,maxind)
  
  probls <- NULL 
  
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
      } else {
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
          if ((fmxind > 2 & fmxind < length(abscaug))){
            zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
          } else if (fmxind <= 2) {
            while (fmxind <= 2) {
              abscaug <- sort(unique(c(seq(abscaug[1],abscaug[2],length.out = 4),abscaug)))
              f <- sapply(abscaug,fbin)
              fmxind <- which(f==max(f))
            }
            zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
          } else if (fmxind==length(abscaug)) {
            while (fmxind==length(abscaug)) {
              abscaug <- sort(unique(c(seq(tail(abscaug,2)[1],tail(abscaug,1),length.out = 3),abscaug)))
              f <- sapply(abscaug,fbin)
              fmxind <- which(f==max(f))
            }
            zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
          }
        } else {
          if (fmxind > 1 & fmxind < (length(abscaug) - 1)) {
            zabsc <- abscaug[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
          } else if (fmxind == 1) {
            while (fmxind==1) {
              abscaug <- sort(unique(c(seq(abscaug[1],abscaug[2],length.out = 3),abscaug)))
              f <- sapply(abscaug,fbin)
              fmxind <- which(f==max(f))
            }
            zabsc <- abscaug[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
          } else if (fmxind > (length(abscaug) - 2)) {
            while (fmxind > (length(abscaug) - 2)) {
              abscaug <- sort(unique(c(seq(tail(abscaug,2)[1],tail(abscaug,1),length.out = 4),abscaug)))
              f <- sapply(abscaug,fbin)
              fmxind <- which(f==max(f))
            }
            zabsc <- abscaug[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
          }
        }
        zval <- intsct2(zabsc,x,n)
        
        abscaug <- sort(unique(c(abscaug,zval)))
        f <- sapply(abscaug,fbin)
        lvec <- fliksum(abscaug,f,zval,maxind)
      }
    }
  }
  meanprobls <- c(meanprobls,mean(probls))
}

xdistls <- sapply(xls, function(x) dbinom(x,n,p0))

meanprobls %*% xdistls
