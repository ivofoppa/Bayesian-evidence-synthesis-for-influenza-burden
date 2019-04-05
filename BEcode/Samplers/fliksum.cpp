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
NumericVector fliksum(DataFrame df, IntegerVector & parms, std::string dist, double zval, int maxind, int zind) {
  double f0,f1,p0,p1,a ;
  NumericVector absc = df[0], f = df[1], lvec ;
  
  for (int k = 1 ; k < absc.size() ; k++) {
    
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


  