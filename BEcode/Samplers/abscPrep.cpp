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
  int maxind = 0, fmxind = 0 ;
  NumericVector f = NumericVector::create(1) ;
  double arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
  double minp = 0, pintv = 0, abscinc = 0, newabsc = 0 ;
  
  absc.sort() ;
  
  f[0] = R::dbinom(x,n,absc[0],true) ;
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
        newabsc = minp + abscinc ;
        absc.insert(1, newabsc) ;
        f.insert(1,R::dbinom(x,n,newabsc,true)) ;
        fmxind = which_max(f) ;
      }
    } else if (fmxind>=(absc.size() - 3)) {
      while (fmxind>=(absc.size() - 3)) {
        pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
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
        newabsc = absc[0] + abscinc ;
        absc.insert(1, newabsc) ;
        f.insert(1,R::dbinom(x,n,newabsc,true)) ;
        fmxind = which_max(f) ;
      }
    } else if (fmxind > (absc.size() - 4)) {
      while (fmxind > (absc.size() - 4)) {
        pintv =  absc[absc.size() - 1] - absc[absc.size() - 2];
        abscinc = pintv/2 ;
        newabsc = absc[absc.size() - 2] + abscinc ;
        absc.insert(absc.size() - 1, newabsc) ;
        f.insert(absc.size() - 1,R::dbinom(x,n,newabsc,true)) ;
        fmxind = which_max(f) ;
      }
    }
  }
  return absc ;
} 
