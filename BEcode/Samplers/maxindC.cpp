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
