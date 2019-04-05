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