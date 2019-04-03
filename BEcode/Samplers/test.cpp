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
NumericVector zabsc2(NumericVector & absc, int x, int n ){
  Function f("zabsc") ;
  NumericVector zabsc = f(absc,x,n) ;
  return zabsc ;
} 

/***R
n <- 100
x <- 1

absc <- c(0,0.02,.05,.1,.15,.21,.25,.4,.5,.6,.7,.8,.9,1)
  zabsc2(absc,x,n) 
*/
