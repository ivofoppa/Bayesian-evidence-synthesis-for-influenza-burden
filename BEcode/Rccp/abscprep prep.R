
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




for (int i=0 ; i < out.length() ; i++) {
  p1 = absc[i] ;
  p2 = absc[i + 1] ;
  p3 = absc[i + 2] ;
  p4 = absc[i + 3] ;
  
  f1 = R::dbinom(x,n,p1,true) ;
  f2 = R::dbinom(x,n,p2,true) ;
  f3 = R::dbinom(x,n,p3,true) ;
  f4 = R::dbinom(x,n,p4,true) ;
  
  a = (f2 - f1)/(p2 - p1) ;
  b = (f4 - f3)/(p4 - p3) ;
  
  p = (f3 - f2 - p3*b + p2*a)/(a - b) ;
  out[i] = p ;
}
return out ;
}

absc <- c(1:10,11:1)/12.

f <- sapply(absc,fbin)
