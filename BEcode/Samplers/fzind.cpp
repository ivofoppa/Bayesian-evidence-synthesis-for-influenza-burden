#include <Rcpp.h>
#include <math.h>
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