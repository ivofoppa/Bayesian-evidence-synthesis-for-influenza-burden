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
// [[Rcpp::export]]
NumericVector pval(NumericVector & absc, int x, int n){
  double a, b ;
  NumericVector zvec(absc.length() - 3), f(absc.length()), pvec(1) ;
  zvec = intsct2(absc,x,n) ;
  
  for (int i = 1 ; i < f.length() ; i++){
    f[i] = R::dbinom(x,n,absc[i], true) ;
  }
  f[0] = -1e+10 ;
  f[f.length() - 1] = -1e+10 ;
  
  for (int i = 0 ; i < (absc.length() - 4) ; i++){
    a = (f[i + 1] - f[i])/(absc[i + 1] - absc[i]) ;
    
    b = (f[i + 3] - f[i + 2])/(absc[i + 3] - absc[i + 2]) ;
    
    pvec.push_back(std::exp(f[i + 1])/a * (std::exp(a * (zvec[i] - absc[i + 1]) - 1))) ;
    pvec.push_back(std::exp(f[i + 2])/b * (1 - std::exp(b * (zvec[i] - absc[i + 2])))) ;
  }
  
  pvec[0] = std::exp(f[1]) * (absc[1] - absc[0])/2 ;
  pvec.push_back(std::exp(f[f.length() - 2]) * (absc[absc.length() - 1] - absc[absc.length() - 2])/2) ;
  
  return pvec ;
}

// [[Rcpp::export]]
int whichP(NumericVector & cumpvec, double x){
  int sumcnt=0 ;
  IntegerVector indvec(sumcnt);
  
  for ( int i=0 ; i < cumpvec.length() ; i++ ){
    if ( x >= cumpvec[i]) {
      indvec.push_back(i) ;
      sumcnt += 1 ;
    }
  }
return indvec.std::end() ;
}

/*** R
absc <- c(1e-15,0.1,.2,.3,.5,.7,1 - 1e-15)
pval(absc,20,100)
  */
