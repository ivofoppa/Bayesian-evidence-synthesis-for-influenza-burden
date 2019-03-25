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
int abscPrep(NumericVector absc, NumericVector f, int x, int n ){
  int maxind = -99 ;
  int fmxind = which_max(f) ;
  double arg1 = std::max(absc[fmxind] - 1e-9,0.) ;
  double arg2 = std::min(absc[fmxind] + 1e-9,1.) ;
  
  if (f[fmxind] > R::dbinom(x,n,arg1,true) && f[fmxind] > R::dbinom(x,n,arg2,true)) {
    maxind = 0 ;
  } else if (f[fmxind] < R::dbinom(x,n,arg2,true)) {
    maxind = 1 ;
  } else if (R::dbinom(x,n,arg1,true) > f[fmxind]) {
    maxind = 2 ;
    }
    return maxind ;
}

/*** R
absc <- c(1e-40,.05,.1,.15,.21,.25,.3,.4,.5,.6,.7,.8,.9,1 - 1e-40)

f <- sapply(absc,function(p) dbinom(30,100,p,TRUE))

abscPrep(absc,f,30,100)
*/

