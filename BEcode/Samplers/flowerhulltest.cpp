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
double flowerhull(double p, DataFrame df, IntegerVector & parms, std::string dist, double zval) {
  double f0,f1,lh ;
  int selind = 0, selind2 = 0 ;
  NumericVector absc = df[0], f = df[1] ;
  
  absc.push_back(zval) ;
  absc.sort() ;

  for (int i = 0 ; i < absc.size() ; i++ ) {
    if (absc[i] < p) {
      selind = i ;
    }
    if (absc[i] <= zval) {
      selind2 = i ;
    }
  }
  
  if (dist=="binom") {  
    int x = parms[0] ;
    int n = parms[1] ;
    f.insert(selind2,R::dbinom(x,n,zval,true)) ;
  } else if (dist=="Poi") {  
    int x = parms[0] ;
    f.insert(selind2,R::dpois(x,zval,true)) ;
   }
  
  double xcap = p - absc[selind] ;
  f0 = f[selind] ;
  f1 = f[selind + 1] ;
  
  if (!std::isinf(f[selind])) {
    lh = exp(f0 + (f1 - f0)/(absc[selind + 1] - absc[selind]) * xcap) ;
  } else {
    lh = 0 ;
  }
  df = DataFrame::create(_["absc"] = absc,_["f"] = f) ;
  
  return lh ;
}
