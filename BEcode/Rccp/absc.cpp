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
int abscPrep(NumericVector absc, NumericVector f ){
  int maxind ;
  NumericVector abscout(2) ;
  
  int fmxind = which_max(absc) ;
  double arg1 = std::max(absc[fmxind] - 1e-5,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-5,1.) ;
  
  if (f[fmxind] > Rcpp::fbin(arg1) && f[fmxind] > Rcpp::fbin(arg2)) {
    maxind = 0 ;
  } else if (f[fmxind] < Rcpp::fbin(arg2))) {
    maxind = 1 ;
  } else if (Rcpp::fbin(arg1) > f[fmxind]) {
    maxind = 2 ;
  }
  return maxind ;
}

/*** R
absc <- c(1:10,11:1)/12.

f <- sapply(absc,fbin)

abscPrep(absc,f)
*/

