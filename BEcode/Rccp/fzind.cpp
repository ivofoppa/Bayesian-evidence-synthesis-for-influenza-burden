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

double lik(double pr, NumericVector & parms, int parmind, DataFrame & data) {
  NumericVector x = data[0], cases=data[1], controls=data[2], lvec ;
  double lv ;
  if (parmind==1) {
    for (int k = 0 ; k < cases.size() ; k++) {
      lv = cases[k]*(pr + parms[1] * x[k] - log(1 + exp(pr + parms[1] * x[k]))) + controls[k]*(- log(1 + exp(pr + parms[1] * x[k])))  ;
      lvec.push_back(lv) ;
    }
  } else if (parmind==2) {
    for (int k = 0 ; k < cases.size() ; k++) {
      lv = cases[k]*(parms[0] + pr * x[k] - log(1 + exp(parms[0] + pr * x[k]))) + controls[k]*(- log(1 + exp(parms[0] + pr * x[k])))  ;
      lvec.push_back(lv) ;
    }
  }
  return sum(lvec) ;
}

// [[Rcpp::export]]
int fzind(DataFrame df, NumericVector & parms, int parmind, DataFrame & data, double zval) {
  double arg ;
  int zind = 0 ;

  arg = zval - 1e-9 ;
  
  if (lik(arg,parms, parmind, data) > lik(zval,parms, parmind, data)) {
    zind = 1 ;
  } else {
    zind = 0 ;
  }
  
  return zind ;
}
