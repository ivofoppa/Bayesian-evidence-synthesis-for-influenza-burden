#include <Rcpp.h>
#include <cmath>
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
double flowerhull(double p, DataFrame & df) {
  int selind = 0 ;
  NumericVector absc = df[0], l = df[1], abscpart ;
  double l0,l1,p0,p1,lh,xcap = 0, pmax ;
  
  abscpart = absc[absc < p] ;
  pmax = max(abscpart) ;
  
  selind = abscpart.length() - 1 ;
  
  xcap = p - pmax ;
  l0 = l[selind] ;
  l1 = l[selind + 1] ;
  
  p0 = absc[selind] ;
  p1 = absc[selind + 1] ;
  
  lh = l0 + (l1 - l0) * xcap/(p1 - p0) ;
  return lh ;
}
