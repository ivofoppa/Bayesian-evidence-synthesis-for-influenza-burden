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
NumericVector intsct2(NumericVector absc, int x, int n){
  double p, f1, f2, f3, f4, p1,p2,p3,p4,a, b ;
  NumericVector out(absc.length() - 3) ;
  for (int i=0 ; i < (absc.length() - 3) ; i++) {
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
    p = (f3 - f1 + p1*a - p3*b)/(a - b) ;
    out[i] = p ;
  }
  return out ;
}

// [[Rcpp::export]]
NumericVector pval(NumericVector absc, int x, int n){
  double f1, f2, f3, f4, a, b ;
  NumericVector intsct(absc.length() - 3), f(absc.length()), pvec(absc.length() - 1) ;
   intsct = intsct2(absc,20,100) ;
   
   for (int i = 0 ; i < (absc.length() - 2) ; i++){
     f[i + 1] = R::dbinom(x,n,absc[i + 1], false) ;
   }
   f[0] = 0 ;
   f[absc.length() - 1] = 0 ;
   
   for (int i = 0 ; i < (absc.length() - 2) ; i++){
     f1 = f[i] ;
     f2 = f[i + 1] ;
     a = (f2 - f1)/(absc[i + 1] - absc[i]) ;
     
     f3 = f[i + 2] ;
     f4 = f[i + 3] ;
     b = (f4 - f3)/(absc[i + 3] - absc[i + 2]) ;
     
     pvec[i + 1] = f1/a * (std::exp(intsct[i] - absc[i + 1]*a) - 1) + 
       f2/b * (1 - std::exp(intsct[i + 1] - absc[i + 2]*b));
   }
   
   pvec[0] = f[1] * (absc[1] - absc[0]) ;
   pvec[absc.length()] = f[absc.length() - 2] * (absc[absc.length() - 1] - absc[absc.length() - 2]) ;
     
  return pvec ;
}

/*** R
absc <- c(1e-10,.1,.15,.2,.25,.4,.5,1.- 1e-10)
pval(absc,20,100)
*/
