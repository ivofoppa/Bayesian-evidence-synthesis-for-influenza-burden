#include <Rcpp.h>
using namespace Rcpp;

DataFrame abscPrep(NumericVector & absc, IntegerVector & parms, std::string & dist ){
  int maxind = 0, fmxind = 0 ;
  NumericVector f ;
  double arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
  double minp = 0, pintv = 0, abscinc = 0, newabsc = 0 ;
  
  absc.sort() ;
  absc = absc[(absc >=0) & (absc <= 1)] ;
  
  if (absc[0] != 0) {
    absc.push_front(0) ;
  }
  if (absc[absc.size() - 1] != 1) {
    absc.push_back(1) ;
  }
  
  if (dist=="binom") {  
    int x = parms[0] ;
    int n = parms[1] ;
    for (int i=0; i < absc.size() ; i++) {
      f.push_back(R::dbinom(x,n,absc[i],true)) ;
    }
    
    fmxind = which_max(f) ;
    
    if (f[fmxind] > R::dbinom(x,n,arg1,true) && f[fmxind] > R::dbinom(x,n,arg2,true)) {
      maxind = 0 ;
    } else if (f[fmxind] < R::dbinom(x,n,arg2,true)) {
      maxind = 1 ;
    } else if (R::dbinom(x,n,arg1,true) > f[fmxind]) {
      maxind = 2 ;
    }
    
    //* calculating intersection x for middle section (containing max); using Rcpp function intsct2
    if  (maxind==2 ) {
      if (fmxind < 3) {
        while (fmxind < 3 ) {
          minp = absc[0];
          pintv = absc[1] - minp ;
          abscinc = pintv/2 ;
          newabsc = minp + abscinc ;
          absc.insert(1, newabsc) ;
          f.insert(1,R::dbinom(x,n,newabsc,true)) ;
          fmxind = which_max(f) ;
        }
      } else if (fmxind>=(absc.size() - 3)) {
        while (fmxind>=(absc.size() - 3)) {
          pintv =  1 - absc[absc.size() - 1];
          abscinc = pintv/2 ;
          newabsc = absc[absc.size() - 2] + abscinc ;
          absc.insert(absc.size() - 1, newabsc) ;
          f.insert(absc.size() - 1, R::dbinom(x,n,newabsc,true)) ;
          fmxind = which_max(f) ;
        }
      } 
    } else {
      
      if (fmxind < 2) {
        while (fmxind < 2) {
          minp = absc[0] ;
          pintv = absc[1] - minp ;
          abscinc = pintv/2 ;
          newabsc = minp + abscinc ;
          absc.insert(1, newabsc) ;
          f.insert(1,R::dbinom(x,n,newabsc,true)) ;
          fmxind = which_max(f) ;
        }
      } else if (fmxind > (absc.size() - 4)) {
        while (fmxind > (absc.size() - 4) ) {
          pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
          abscinc = pintv/2 ;
          newabsc = absc[absc.size() - 2] + abscinc ;
          absc.insert(absc.size() - 1, newabsc) ;
          f.insert(absc.size() - 1,R::dbinom(x,n,newabsc,true)) ;
          fmxind = which_max(f) ;
        }
      }
    }
  } else if (dist=="Poi") {  
    
    int x = parms[0] ;
    
    for (int i = 0 ; i < absc.size() ; i++) {
      absc[i] = absc[i] * x ;
    }
    
    for (int i=0; i < absc.size() ; i++) {
      f.push_back(R::dpois(x,absc[i],true)) ;
    }
    
    fmxind = which_max(f) ;
    double fabscmin = min(f) ;
    
    if (f[fmxind] > R::dpois(x,arg1,true) && f[fmxind] > R::dpois(x,arg2,true)) {
      maxind = 0 ;
    } else if (f[fmxind] < R::dpois(x,arg2,true)) {
      maxind = 1 ;
    } else if (R::dpois(x,arg1,true) > f[fmxind]) {
      maxind = 2 ;
    }
    
    //* calculating intersection x for middle section (containing max); using Rcpp function intsct2
    if  (maxind==2 ) {
      if (fmxind < 3) {
        while (fmxind < 3) {
          minp = absc[0];
          pintv = absc[1] - minp ;
          abscinc = pintv/2 ;
          newabsc = minp + abscinc ;
          absc.insert(1, newabsc) ;
          f.insert(1,R::dpois(x,newabsc,true)) ;
          fmxind = which_max(f) ;
        }
      } 
      if (fmxind>=(absc.size() - 3) || fabscmin > -50) {
        while (fmxind>=(absc.size() - 3) || fabscmin > -50) {
          pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
          abscinc = pintv * 2 ;
          newabsc = absc[absc.size() - 1] + abscinc ;
          absc.push_back(newabsc) ;
          f.push_back(R::dpois(x,newabsc,true)) ;
          fmxind = which_max(f) ;
          fabscmin = f[f.size() - 1] ;
        }
      }
    } else {
      
      if (fmxind < 2) {
        while (fmxind < 2) {
          minp = absc[0] ;
          pintv = absc[1] - minp ;
          abscinc = pintv/2 ;
          newabsc = absc[0] + abscinc ;
          absc.insert(1, newabsc) ;
          f.insert(1,R::dpois(x,newabsc,true)) ;
          fmxind = which_max(f) ;
        }
      } else if (fmxind > (absc.size() - 4) || fabscmin > -50) {
        while (fmxind > (absc.size() - 4) || fabscmin > -50) {
          pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
          abscinc = pintv * 2 ;
          newabsc = absc[absc.size() - 1] + abscinc ;
          absc.push_back(newabsc) ;
          f.push_back(R::dpois(x,newabsc,true)) ;
          fmxind = which_max(f) ;
          fabscmin = f[f.size() - 1] ;
        }
      }
    }
  }
  NumericVector absc2 ;
  for (int i=0 ; i < absc.size() ; i++){
    if (f[i] > -50 || i==0 || i==(absc.size() - 1)) {
      absc2.push_back(absc[i]);
    } 
  }
  NumericVector f2 ;
  for (int i=0 ; i < absc.size() ; i++){
    if (f[i] > -50 || i==0 || i==(absc.size() - 1)) {
      f2.push_back(f[i]);
    } 
  }
  DataFrame df = DataFrame::create( Named("absc")=absc2, Named("f")=f2 ) ;
  return df;
}

int maxindC(DataFrame & df, IntegerVector & parms, std::string & dist ){
  int maxind = 0, fmxind = 0 ;
  double arg1, arg2 ;
  NumericVector absc = df[0], f = df[1] ;
  if (dist=="binom") {  
    int x = parms[0] ;
    int n = parms[1] ;
    
    fmxind = which_max(f) ;
    
    arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
    arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
    
    if (f[fmxind] > R::dbinom(x,n,arg1,true) && f[fmxind] > R::dbinom(x,n,arg2,true)) {
      maxind = 0 ;
    } else if (f[fmxind] < R::dbinom(x,n,arg2,true)) {
      maxind = 1 ;
    } else if (R::dbinom(x,n,arg1,true) > f[fmxind]) {
      maxind = 2 ;
    }
  } else if (dist=="Poi") {  
    
    int x = parms[0] ;
    
    fmxind = which_max(f) ;
    
    arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
    arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
    
    if (f[fmxind] > R::dpois(x,arg1,true) && f[fmxind] > R::dpois(x,arg2,true)) {
      maxind = 0 ;
    } else if (f[fmxind] < R::dpois(x,arg2,true)) {
      maxind = 1 ;
    } else if (R::dpois(x,arg1,true) > f[fmxind]) {
      maxind = 2 ;
    }
  }
  return maxind ;
} 

NumericVector zabscPrep(DataFrame & df, IntegerVector & parms, std::string & dist, int maxind ){
  
  NumericVector zabsc(4), absc = df[0], f = df[1] ;
  int fmxind = which_max(f) ;
  
  if (maxind == 2) {
    for (int i=0 ; i < 4 ; i++) {
      zabsc[i] = absc[(fmxind - 2 + i)] ;
    }
  } else {
    for (int i=0 ; i < 4 ; i++) {
      zabsc[i] = absc[(fmxind - 1 + i)] ;
    } 
  }
  return zabsc ;
}

double intsct(NumericVector zabsc, IntegerVector & parms, std::string & dist){
  double p = 1,p1,p2,p3,p4,a = 1,b = 0,f1 = 1,f2 = 1,f3 = 1,f4 = 1 ;
  NumericVector out(zabsc.length() - 3) ;
  
  p1 = zabsc[0] ;
  p2 = zabsc[1] ;
  p3 = zabsc[2] ;
  p4 = zabsc[3] ;
  
  if (dist=="binom") {
    int x = parms[0], n = parms[1] ;
    
    f1 = R::dbinom(x,n,p1,true) ;
    f2 = R::dbinom(x,n,p2,true) ;
    f3 = R::dbinom(x,n,p3,true) ;
    f4 = R::dbinom(x,n,p4,true) ;
    
    a = (f2 - f1)/(p2 - p1) ;
    b = (f4 - f3)/(p4 - p3) ;
    
  } else if (dist=="Poi") {
    int x = parms[0] ;
    
    f1 = R::dpois(x,p1,true) ;
    f2 = R::dpois(x,p2,true) ;
    f3 = R::dpois(x,p3,true) ; 
    f4 = R::dpois(x,p4,true) ;
    
    a = (f2 - f1)/(p2 - p1) ;
    b = (f4 - f3)/(p4 - p3) ;
  }
  p = (f3 - f2 - p3*b + p2*a)/(a - b) ;
  
  return p ;
}

DataFrame fabscaug(DataFrame df, IntegerVector & parms, std::string dist, double zval) {
  int selind = 0 ;
  NumericVector absc = df[0], f = df[1] ;
  
  absc.push_back(zval) ;
  absc = unique(absc) ;
  absc.sort() ;
  
  for (int i = 0 ; i < absc.size() ; i++ ) {
    if (absc[i] <= zval) {
      selind = i ;
    }
  }
  
  if (dist=="binom") {  
    int x = parms[0] ;
    int n = parms[1] ;
    f.insert(selind,R::dbinom(x,n,zval,true)) ;
  } else if (dist=="Poi") {  
    int x = parms[0] ;
    f.insert(selind,R::dpois(x,zval,true)) ;
  }
  return DataFrame::create(_["absc"] = absc,_["f"] = f) ;
}

int fzind(DataFrame df, IntegerVector & parms, std::string dist, double zval) {
  double arg ;
  int zind = 0 ;
  NumericVector absc = df[0], f = df[1] ;
  
  arg = std::max(zval - 1e-9,zval/2) ;
  
  if (dist=="binom") {  
    int x = parms[0] ;
    int n = parms[1] ;
    
    if (R::dbinom(x,n,arg,true) > R::dbinom(x,n,zval,true)) {
      zind = 1 ;
    } else {
      zind = 0 ;
    }
  } else if (dist=="Poi") {  
    
    int x = parms[0] ;
    
    if (R::dpois(x,arg,true) > R::dpois(x,zval,true)) {
      zind = 1 ;
    } else {
      zind = 0 ;
    }
  }
  return zind ;
}

NumericVector fliksum(DataFrame df, IntegerVector & parms, std::string dist, double zval, int maxind, int zind) {
  double f0,f1,p0,p1,a ;
  NumericVector absc = df[0], f = df[1], lvec ;
  
  for (int k = 0 ; k < (absc.size() - 1) ; k++) {
    
    if ((absc[k + 1]==zval && zind==0) || (absc[k + 1] < zval )) {
      lvec.push_back(exp(f[k + 1])*(absc[k + 1] - absc[k])) ;
    } else if (absc[k + 1] == zval && zind==1 && maxind!=0){
      f0 = f[k - 1] ;
      f1 = f[k] ;
      p0 = absc[k - 1] ;
      p1 = absc[k] ;
      a = (f1 - f0)/(p1 - p0) ;
      lvec.push_back(exp(f1)*(exp(a*(zval - p1)) - 1)/a) ;
    } else if (absc[k] == zval && zind==0 && maxind!=0) {
      f0 = f[k + 1] ;
      f1 = f[k + 2] ;
      p0 = absc[k + 1] ;
      p1 = absc[k + 2] ;
      a = (f1 - f0)/(p1 - p0) ;
      lvec.push_back(exp(f0)*(1 - exp(a*(zval - p0)))/a) ;
    } else if ((absc[k]==zval && zind==1) || (absc[k + 1]==zval && zind==1 && maxind==0) || 
      (absc[k]==zval && zind==0 && maxind==0) || (absc[k] > zval)) {
      lvec.push_back(exp(f[k])*(absc[k + 1] - absc[k])) ;
    } 
  }
  return lvec ;
}

double fpsample(DataFrame df, IntegerVector & parms, std::string dist, double zval, int maxind, int zind, NumericVector lvec) {
  double f0,f1,p0,p1,a,pran,pbase,ptop,prm,pout ;
  int selind = 0 ;
  
  NumericVector absc = df[0], f = df[1], pvec, pveccum ;
  
  for (int i = 0 ; i < lvec.size() ; i++) {
    pvec.push_back(lvec[i]/sum(lvec)) ;
    if (i > 0) {
      pveccum.push_back(pveccum[i - 1] + pvec[i]) ;
    } else pveccum.push_back(pvec[i]) ;
  }
  
  pveccum.push_front(0) ;
  
  pran = R::runif(0,1) ;
  
  for (int i = 0 ; i < pveccum.size() ; i++ ) {
    if (pveccum[i] <= pran) {
      selind = i ;
    }
  }
  
  pbase = pveccum[selind] ;
  ptop = pveccum[selind + 1] ;
  prm = pran - pbase ;
  
  if (absc[selind + 1] == zval && zind==1 && maxind > 0) {
    f0 = f[selind - 1] ;
    f1 = f[selind] ;
    p0 = absc[selind - 1] ;
    p1 = absc[selind] ;
    a = (f1 - f0)/(p1 - p0) ;
    pout = (log(exp(f1) + a*prm*sum(lvec)) + a*p1 - f1)/a ;
  } else if (absc[selind] == zval && zind==0 && maxind > 0) {
    f0 = f[selind + 1] ;
    f1 = f[selind + 2] ;
    p0 = absc[selind + 1] ;
    p1 = absc[selind + 2] ;
    a = (f1 - f0)/(p1 - p0) ;
    pout = (log(exp(f0 + a*(zval - p0)) + prm*a*sum(lvec)) + a*p0 - f0)/a ;
  } else {
    pout = absc[selind] + (absc[selind + 1] - absc[selind])*prm/(ptop - pbase) ;
  } 
  return pout ;
}

double flowerhull(double p, DataFrame df, IntegerVector & parms) {
  double f0,f1,lh ;
  int selind = 0 ;
  NumericVector absc = df[0], f = df[1] ;
  
  for (int i = 0 ; i < absc.size() ; i++ ) {
    if (absc[i] < p) {
      selind = i ;
    }
  }
  
  double xcap = p - absc[selind] ;
  f0 = f[selind] ;
  f1 = f[selind + 1] ;
  
  if (!std::isinf(f[selind])) {
    lh = exp(f0 + (f1 - f0)/(absc[selind + 1] - absc[selind]) * xcap) ;
  } else {
    lh = 0 ;
  }
  return lh ;
}

double fupperhull(double p, DataFrame df, IntegerVector & parms, std::string dist, double zval, int maxind, int zind) {
  double f0,f1,uh = 0,p0,p1,a ;
  int selind = 0 ;
  NumericVector absc = df[0], f = df[1] ;
  
  for (int i = 0 ; i < absc.size() ; i++ ) {
    if (absc[i] < p) {
      selind = i ;
    }
  }
  
  if ((absc[selind + 1]==zval && zind==0) || (absc[selind + 1] < zval )){
    uh = exp(f[selind + 1]) ;
  } else if (absc[selind + 1] == zval && zind==1 && maxind!=0){
    f0 = f[selind - 1] ;
    f1 = f[selind] ;
    p0 = absc[selind - 1] ;
    p1 = absc[selind] ;
    a = (f1 - f0)/(p1 - p0) ;
    uh = exp(f1 + a*(p - p1)) ;
  } else if (absc[selind] == zval && zind==0 && maxind!=0){
    f0 = f[selind + 1] ;
    f1 = f[selind + 2] ;
    p0 = absc[selind + 1] ;
    p1 = absc[selind + 2] ;
    a = (f1 - f0)/(p1 - p0) ;
    uh = exp(f0 + a*(p - p0)) ;
  } else if ((absc[selind]==zval && zind==1) || (absc[selind + 1]==zval && zind==1 && maxind==0) || 
    (absc[selind]==zval && zind==0 && maxind==0) || (absc[selind] > zval)) {
    uh = exp(f[selind]) ;
  }
  return uh ;
}

// [[Rcpp::export]]
NumericVector samples(NumericVector & absc, IntegerVector & parms, std::string dist, int nsims, double crit) {
  int maxind = -99, zind = -99, x, n ;
  double zval, fval, uran, psample, lh, uh, uratio ;
  NumericVector zabsc, lvec, parmsmplvec, f ;
  DataFrame df ;
  
  while (parmsmplvec.size() < nsims ) {
    df = abscPrep(absc,parms,dist) ;
    maxind = maxindC(df,parms,dist) ;
    zabsc = zabscPrep(df,parms,dist,maxind) ;
    zval = intsct(zabsc,parms, dist) ;
    df = fabscaug(df,parms,dist,zval) ;
    f = df[1] ;
    zind = fzind(df,parms,dist,zval) ;
    lvec = fliksum(df, parms, dist, zval, maxind,zind) ;
    
    psample = fpsample(df, parms, dist, zval, maxind,zind,lvec) ;
    uran = R::runif(0,1) ;
    lh = flowerhull(psample,df,parms) ;
    uh = fupperhull(psample,df,parms,dist,zval,maxind,zind) ;
    if (uran < lh/uh) {
      parmsmplvec.push_back(psample) ;
    } else {
      if (dist=="binom") {
        x = parms[0] ;
        n = parms[1] ;
        fval = R::dbinom(x,n,psample,false) ;
      } else {
        x = parms[0] ;
        fval = R::dpois(x,psample,false) ;
        }
      uratio = fval/uh ;
      if (uran < uratio) {
        parmsmplvec.push_back(psample) ;
      } else if (uratio < crit ) {
        absc.push_back(psample) ;
        absc.sort() ;
      }
    }
  }
  return parmsmplvec ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
parms <- c(10,100)
dist <- "Poi"
dist <- "binom"
absc <- c(.05,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9)

start_time <- Sys.time()
sls <- samples(absc,parms,dist,10000,.75)
end_time <- Sys.time()
end_time - start_time
hist(sls,100,xlim = c(0,.3),freq = F)
lines(ls,yls*100,col = 'red')
*/
