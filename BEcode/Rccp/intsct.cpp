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

double lik(double pr,NumericVector & parms,int parmind,DataFrame & data) {
  NumericVector x = data[0],cases=data[1],controls=data[2],lvec ;
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
double intsct(DataFrame & df,NumericVector & parms,int parmind,DataFrame & data,int maxind){
  
  NumericVector zabsc(4),absc = df[0],f = df[1] ;
  int fmxind = which_max(f) ;
  double p = 1,p1,p2,p3,p4,a = 1,b = 0,f1 = 1,f2 = 1,f3 = 1,f4 = 1 ;
  NumericVector out(zabsc.length() - 3) ;
  
  if (maxind == 2) {
    for (int i=0 ; i < 4 ; i++) {
      zabsc[i] = absc[(fmxind - 2 + i)] ;
    }
  } else {
    for (int i=0 ; i < 4 ; i++) {
      zabsc[i] = absc[(fmxind - 1 + i)] ;
    } 
  }
  
  p1 = zabsc[0] ;
  p2 = zabsc[1] ;
  p3 = zabsc[2] ;
  p4 = zabsc[3] ;
  
  f1 = lik(p1,parms,parmind,data);
  f2 = lik(p2,parms,parmind,data);
  f3 = lik(p3,parms,parmind,data);
  f4 = lik(p4,parms,parmind,data);
  
  a = (f2 - f1)/(p2 - p1) ;
  b = (f4 - f3)/(p4 - p3) ;
  
  p = (f3 - f2 - p3*b + p2*a)/(a - b) ;
  
  return p ;
}

