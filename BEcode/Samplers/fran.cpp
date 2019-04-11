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
double lik(double pr, int parmind, NumericVector & parms, DataFrame & data) {
  NumericVector c = data[0], x=data[1], cnt=data[2] ;
  double lv = 1 ;
  if (parmind==1) {
    for (int i = 0; i < x.length(); i++) {
      lv *= pow(exp(pr + parms[1] * x[i]),c[i])/(1 + exp(pr + parms[1] * x[i])) ;
    }
  } else if (parmind==2) {
    for (int i = 0; i < x.length(); i++) {
      lv *= pow(exp(parms[0] + pr * x[i]),c[i])/(1 + exp(parms[0] + pr * x[i])) ;
    }
  }
  
  return lv ;
} 
