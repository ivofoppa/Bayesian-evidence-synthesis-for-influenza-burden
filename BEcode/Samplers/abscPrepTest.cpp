#include <Rcpp.h>
#include <limits>
#include <vector>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
  //
  //   http://www.rcpp.org/
  //   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
  //
  // [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
NumericVector abscPrep(NumericVector & absc, int x, int n ){
  int maxind = -99, fmxind = 0, minpexp = -5 ;
  NumericVector f(0) ;
  double arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
  double minp = 0, pintv = 0, abscinc = 0 ;
  
  absc.sort() ;
  
  f = (R::dbinom(x,n,absc[0],true)) ;
  for (int i=1; i < absc.size() ; i++) {
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
      while (fmxind < 3) {
        minp = absc[0];
        pintv = absc[1] - minp ;
        abscinc = pintv/2 ;
        absc.insert(1, minp + abscinc) ;
        absc.sort() ;
        f = (R::dbinom(x,n,absc[0],true)) ;
        for (int i=1 ; i < absc.size() ; i++) {
          f.push_back(R::dbinom(x,n,absc[i],true)) ;
        }
        fmxind = which_max(f) ;
      }
    } else if (fmxind>=(absc.size() - 3)) {
      while (fmxind>=(absc.size() - 3)) {
        pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
        abscinc = pintv/2 ;
        absc.insert(absc.size() - 2, absc[absc.size() - 2] + abscinc) ;
        absc.sort() ;
        f = (R::dbinom(x,n,absc[0],true)) ;
        for (int i=1; i < absc.size() ; i++) {
          f.push_back(R::dbinom(x,n,absc[i],true)) ;
        }
        fmxind = which_max(f) ;
      }
    }
  } else {
    
    if (fmxind < 2) {
      while (fmxind < 2) {
        minp = absc[0] ;
        pintv = absc[1] - minp ;
        abscinc = pintv/2 ;
        absc.insert(1, absc[1] + abscinc) ;
        absc.sort() ;
        f = (R::dbinom(x,n,absc[0],true)) ;
        for (int i=1 ; i < absc.size() ; i++) {
          f.push_back(R::dbinom(x,n,absc[i],true)) ;
        }
        fmxind = which_max(f) ;
        minpexp = minpexp*2 ;
      }
    } else if (fmxind > (absc.size() - 4)) {
      while (fmxind > (absc.size() - 4)) {
        pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
        abscinc = pintv/2 ;
        absc.insert(absc.size() - 2, absc[absc.size() - 2] + abscinc) ;
        absc.sort() ;
        f = (R::dbinom(x,n,absc[0],true)) ;
        for (int i=1; i < absc.size() ; i++) {
          f.push_back(R::dbinom(x,n,absc[i],true)) ;
        }
        fmxind = which_max(f) ;
      }
    }
  }
  return absc ;
} 
// [[Rcpp::export]]
NumericVector zabsc(NumericVector & absc0, int x, int n ){
  double arg1 = 0, arg2 = 0 ;
  
  int fmxind = 0, maxind = 0 ;
  NumericVector zabsc(4), absc(0), f(0) ;
  
  absc = abscPrep(absc0,x,n) ;
  f = (R::dbinom(x,n,absc[0],true)) ;
  for (int i=1 ; i < absc.size() ; i++) {
    f.push_back(R::dbinom(x,n,absc[i],true)) ;
  }
  
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

/*** R
absc <- c(0,0.02,.05,.1,.15,.21,.25,.3,.4,.5,.6,.7,.8,.9,1)

x <- 1 ; n <- 100
f <- sapply(absc,function(p) dbinom(x,n,p,TRUE))

zabsc(absc,x,n)
  */

