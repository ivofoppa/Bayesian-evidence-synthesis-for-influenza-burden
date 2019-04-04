#include <Rcpp.h>
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

// [[Rcpp::export]]
NumericVector zabscPrep(NumericVector & absc, IntegerVector & parms, std::string & dist ){
  double arg1 = 0, arg2 = 0 ;
  
  int fmxind = 0, maxind = 0 ;
  NumericVector zabsc(4), f(0) ;
  
  if (dist=="binom"){
    int x = parms[0], n = parms[1] ;
    
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
  } else if (dist=="Poi"){
    int x = parms[0] ;
    
    f = (R::dpois(x,absc[0],true)) ;
    for (int i=1 ; i < absc.size() ; i++) {
      f.push_back(R::dpois(x,absc[i],true)) ;
    }
    
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