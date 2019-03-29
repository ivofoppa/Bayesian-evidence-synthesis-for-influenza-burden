#include <Rcpp.h>
#include <limits>
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
NumericVector intsct2(NumericVector absc, int x, int n){
  double p, f1, f2, f3, f4, p1,p2,p3,p4,a, b ;
  NumericVector out(absc.length() - 3) ;
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