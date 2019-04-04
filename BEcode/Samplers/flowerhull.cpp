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
double flowerhull(double p,NumericVector & absc, IntegerVector & parms, std::string & dist) {
  double f0,f1,lh ;
  int selind = 0 ;
  NumericVector f ;
  
  if (dist=="binom") {
    int x = parms[0], n = parms[1] ;
    
    for (int i=0 ; i < absc.size() ; i++) {
      f.push_back(R::dbinom(x,n,absc[i],true)) ;
    }
  } else if (dist=="Poi") {
    int x = parms[0] ;
    
    for (int i=0 ; i < absc.size() ; i++) {
      f.push_back(R::dpois(x,absc[i],true)) ;
    }
  }  

  for (int i = 0 ; i < absc.size() ; i++ ) {
    if (absc[i] < p) {
      selind = i ;
    }
  }
  
  double xcap = p - absc[selind] ;
  f0 = f[selind] ;
  f1 = f[selind + 1] ;
  
  if (!std::isinf(f[selind])) {
    lh = exp(f0 + (f1 - f0)/(absc[selind + 1] - absc[selind]) * xcap) ;
  } else {
    lh = 0 ;
  }
  return lh ;
}
