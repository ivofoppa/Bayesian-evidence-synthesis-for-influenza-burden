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

double lik1(double pr, int parmind, NumericVector & parms, DataFrame & data) {
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

double lik2(double pr, int parmind, NumericVector & parms, DataFrame & data) {
  NumericVector x = data[0], cases=data[1], pt=data[2], factvec, lvec ;
  double lv ;
  factvec = factorial(cases) ;
  if (parmind==1) {
    for (int k = 0 ; k < cases.size() ; k++) {
      lv =  -exp(pr + parms[1] * x[k] + log(pt[k])) + cases[k] *(pr + parms[1] * x[k] + log(pt[k])) - log(factvec[k]);
      lvec.push_back(lv) ;
    }
  } else if (parmind==2) {
    for (int k = 0 ; k < cases.size() ; k++) {
      lv =  -exp(parms[0] + pr * x[k] + log(pt[k])) + cases[k] * (parms[0] + pr * x[k] + log(pt[k])) - log(factvec[k]);
      lvec.push_back(lv) ;
    }
  }
  return sum(lvec) ;
}

// [[Rcpp::export]]
SEXP callFun(double pr, int parmind, NumericVector & parms, DataFrame & data, Function fun) {

  SEXP lv ; 
  lv = fun(pr, parmind, parms, data) ;
  
  return lv ;
}

/***R
N0 <- 10000000
N1 <- N0
  
Nls <- c(N0,N1)
  
phi <- 2;

lam0 <- 1000
lamls <- c(lam0,lam0*phi)
  
casels <- rpois(2,lamls)
  
  data <- data.frame(exp = c(0,1), cases = casels, pt = Nls)
  
  parmind <- 2

parms <- c(1e-4,0)
  
callFun(.0001, parmind, parms, data,lik2)
*/
