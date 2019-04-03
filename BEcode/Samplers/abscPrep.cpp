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
NumericVector abscPrep(NumericVector & absc, IntegerVector & parms, std::string & dist ){
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
    return absc2 ;
} 
