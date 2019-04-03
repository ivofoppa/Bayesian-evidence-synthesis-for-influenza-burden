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
double fupperhull(double p,NumericVector & absc, int x, int n ) {
  double f0,f1,uh ;
  int selind = 0, zind = 0 ;
  NumericVector f(0) ;
  
  f = (R::dbinom(x,n,absc[0],true)) ;
  for (int i=1 ; i < absc.size() ; i++) {
    f.push_back(R::dbinom(x,n,absc[i],true)) ;
  }
  
  for (int i = 0 ; i < absc.size() ; i++ ) {
    if (absc[i] < p) {
      selind = i ;
    }
  }
  
  double xcap = p - absc[selind] ;
  if (fdist(max(zval - 1e-5,zval/2)) > fdist(zval) , 1, 0)
    if ((absc[selind + 1]==zval & zind==0) | (absc[selind + 1]<zval )){
      fval <- exp(f[selind + 1])
    } else if (absc[selind + 1] == zval & zind==1 & maxind!=0){
      f0 <- f[selind - 1]
      f1 <- f[selind]
      p0 <- absc[selind - 1]
      p1 <- absc[selind]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f1 + a*(p - p1))
    } else if (absc[selind] == zval & zind==0 & maxind!=0){
      f0 <- f[selind + 1]
      f1 <- f[selind + 2]
      p0 <- absc[selind + 1]
      p1 <- absc[selind + 2]
      a <- (f1 - f0)/(p1 - p0)
      fval <- exp(f0 + a*(p - p0))
    } else if ((absc[selind]==zval & zind==1) | (absc[selind + 1]==zval & zind==1 & maxind==0) | 
      (absc[selind]==zval & zind==0 & maxind==0) | (absc[selind] > zval)) {
      fval <- exp(f[selind])
    } 
    fval
  }
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//