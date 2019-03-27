setwd('C:/Users/vor1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/BEcode/Samplers')
source('binom derivative free 5.R')

n <- 100
x <- 20

nsim <- 10000

absc <- c(0,0.01,.05,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9,.99,1)


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

### calculating intersection x for middle section (containing max); using Rcpp function intsct2
if (maxind==2) {
  if ((fmxind > 3 & fmxind < (length(absc) - 1))){
    zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
  } else if (fmxind <= 3) {
    while (fmxind <= 3) {
      minp <- absc[1]
      absc <- sort(unique(c(seq(minp,absc[2],length.out = 3),absc)))
      f <- sapply(absc,fbin)
      fmxind <- which(f==max(f))
    }
    zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
  } else if (fmxind >= (length(absc) - 1)) {
    while (fmxind >= (length(absc) - 1)) {
      maxp <- 1
      absc <- sort(unique(c(seq(tail(absc,2)[1],tail(absc,1),length.out = 3),absc)))
      f <- sapply(absc,fbin)
      fmxind <- which(f==max(f))
    }
    zabsc <- absc[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
  }
} else {
  if (fmxind > 2 & fmxind < (length(absc) - 2)) {
    zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
  } else if (fmxind <= 2) {
    while (fmxind <= 2) {
      absc <- sort(unique(c(seq(absc[1],absc[2],length.out = 3),absc)))
      f <- sapply(absc,fbin)
      fmxind <- which(f==max(f))
    }
    zabsc <- absc[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
  } else if (fmxind > (length(absc) - 3)) {
    while (fmxind > (length(absc) - 3)) {
      absc <- sort(unique(c(seq(tail(absc,2)[1],tail(absc,1),length.out = 3),absc)))
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

dx <- 0.0001
ls <- seq(0.0001,.999999999,dx)
# ls <- seq(0.2,.7,dx)

uhls <- sapply(ls, function(p) fupperhull(p,abscaug,f,zval,maxind))/(sum(exp(yls))*dx)
# uhls <- sapply(ls, function(p) fupperhull2(p,abscaug,f,zval,maxind))
yls <- sapply(ls,fbin)

lhls <- sapply(ls, function(p) flowerhull(p,abscaug,f))/(sum(exp(yls))*dx)

muhls <- max(uhls)
setwd('C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/BEwriteup/Graphics')

file.pdf <- paste0('adaptation_sr',1,'.pdf')
pdf(file.pdf,paper='USr') 
plot(ls,exp(yls)/(sum(exp(yls))*dx),type = 'l',xlim=c(0,.4),ylim = c(0,muhls),xlab = 'p',ylab = '')
# plot(ls,exp(yls),type = 'l',xlim=c(.2,.8),ylim = c(0,.2))
lines(ls,lhls,col = 'blue')
lines(ls,uhls,col = 'red')
dev.off()

adaptn <- 2
adlist <- c(1:10,seq(20,150,10),300,500,1000)

crit <- .75

probls <- NULL
while (adaptn <= 1000) {
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
    } else if (aratio < crit) {
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
        if ((fmxind > 3 & fmxind < (length(abscaug) - 1))){
          zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
        } else if (fmxind <= 3) {
          while (fmxind <= 3) {
            minp <- abscaug[1]
            abscaug <- sort(unique(c(seq(minp,abscaug[2],length.out = 3),abscaug)))
            f <- sapply(abscaug,fbin)
            fmxind <- which(f==max(f))
          }
          zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
        } else if (fmxind >= (length(abscaug) - 1)) {
          while (fmxind >= (length(abscaug) - 1)) {
            maxp <- 1
            abscaug <- sort(unique(c(seq(tail(abscaug,2)[1],tail(abscaug,1),length.out = 3),abscaug)))
            f <- sapply(abscaug,fbin)
            fmxind <- which(f==max(f))
          }
          zabsc <- abscaug[c(fmxind - 2,fmxind - 1,fmxind,fmxind + 1)]
        }
      } else {
        if (fmxind > 2 & fmxind < (length(abscaug) - 2)) {
          zabsc <- abscaug[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
        } else if (fmxind <= 2) {
          while (fmxind <= 2) {
            abscaug <- sort(unique(c(seq(abscaug[1],abscaug[2],length.out = 3),abscaug)))
            f <- sapply(abscaug,fbin)
            fmxind <- which(f==max(f))
          }
          zabsc <- abscaug[c(fmxind - 1,fmxind,fmxind + 1,fmxind + 2)]
        } else if (fmxind > (length(abscaug) - 3)) {
          while (fmxind > (length(abscaug) - 3)) {
            abscaug <- sort(unique(c(seq(tail(abscaug,2)[1],tail(abscaug,1),length.out = 3),abscaug)))
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
      
      if (adaptn%in%adlist) {
        zval <- intsct2(zabsc,x,n)
        
        abscaug <- sort(unique(c(abscaug,zval)))
        f <- sapply(abscaug,fbin)
        lvec <- fliksum(abscaug,f,zval,maxind)
        
        uhls <- sapply(ls, function(p) fupperhull(p,abscaug,f,zval,maxind))/(sum(exp(yls))*dx)
        # uhls <- sapply(ls, function(p) fupperhull2(p,abscaug,f,zval,maxind))
        yls <- sapply(ls,fbin)
        
        lhls <- sapply(ls, function(p) flowerhull(p,abscaug,f))/(sum(exp(yls))*dx)
        
        setwd('C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/BEwriteup/Graphics')
        
        file.pdf <- paste0('adaptation_sr',adaptn,'.pdf')
        pdf(file.pdf,paper='USr') 
        plot(ls,exp(yls)/(sum(exp(yls))*dx),type = 'l',xlim=c(0,.4),ylim = c(0,muhls),xlab = 'p',ylab = '')
        # plot(ls,exp(yls),type = 'l',xlim=c(.2,.8),ylim = c(0,.2))
        lines(ls,lhls,col = 'blue')
        lines(ls,uhls,col = 'red')
        dev.off()
      }
      adaptn <- adaptn + 1
    }
  }
}

