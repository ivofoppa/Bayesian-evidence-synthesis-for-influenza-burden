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
DataFrame fabscaug(DataFrame df, IntegerVector & parms, std::string dist, double zval) {
  int selind = 0 ;
  NumericVector absc = df[0], f = df[1] ;
  
  absc.push_back(zval) ;
  absc = unique(absc) ;
  absc.sort() ;

  for (int i = 0 ; i < absc.size() ; i++ ) {
    if (absc[i] <= zval) {
      selind = i ;
    }
  }
  
  if (dist=="binom") {  
    int x = parms[0] ;
    int n = parms[1] ;
    f.insert(selind,R::dbinom(x,n,zval,true)) ;
  } else if (dist=="Poi") {  
    int x = parms[0] ;
    f.insert(selind,R::dpois(x,zval,true)) ;
   }
  return DataFrame::create(_["absc"] = absc,_["f"] = f) ;
}
